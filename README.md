# 双端 FASTQ 质控与覆盖核查 API

FastAPI 实现的成对 FASTQ 质量过滤服务：上传 `R1.fastq.gz` / `R2.fastq.gz`
与 Phred+33 质量阈值，后台流式处理，支持查询状态、停止任务、下载过滤后的
成对 gzip 结果。

在质控结果之上还提供覆盖核查：提交一个成功质控任务号、参考序列 FASTA、
片段长度范围与最低深度，后台流式重读已发布的 R1/R2，返回需要补测的低覆盖
区间和三类读段对数。

## 快速启动

```bash
docker compose up --build
```

启动后：

* Swagger 文档：<http://localhost:8000/docs>
* OpenAPI JSON：<http://localhost:8000/openapi.json>
* 健康检查：<http://localhost:8000/health>

数据目录 `/data`（Compose 卷 `qc-data`）保存上传文件和结果；容器重启后
内存中的任务记录会失效，遗留目录会被自动清理。

不装 Docker 时也可本地运行：

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
QC_DATA_DIR=./data uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 接口一览

| 方法 | 路径 | 说明 |
|------|------|------|
| `POST` | `/jobs` | 上传 `r1`、`r2`（multipart `File`）和 `min_quality`，返回 `202` 与任务号 |
| `GET` | `/jobs/{job_id}` | 查询状态、进度、保留/丢弃对数、错误信息、下载链接 |
| `GET` | `/jobs` | 列出所有任务（新到旧） |
| `POST` | `/jobs/{job_id}/stop` | 请求停止；已完成/失败的任务返回 `409` |
| `GET` | `/jobs/{job_id}/download/r1` | 下载过滤后的 R1（仅成功任务） |
| `GET` | `/jobs/{job_id}/download/r2` | 下载过滤后的 R2（仅成功任务） |
| `POST` | `/coverage-jobs` | 提交覆盖核查：`qc_job_id`、`reference`（FASTA）、`min_fragment`、`max_fragment`、`min_depth`，返回 `202` 与任务号 |
| `POST` | `/jobs/{job_id}/coverage` | 同上，质控任务号走路径，其余字段走表单 |
| `GET` | `/coverage-jobs` | 列出覆盖任务（新到旧） |
| `GET` | `/coverage-jobs/{job_id}` | 查询覆盖任务状态与报告（成功后含低覆盖区间） |
| `POST` | `/coverage-jobs/{job_id}/stop` | 请求停止覆盖任务 |
| `GET` | `/coverage-jobs/{job_id}/report` | 下载已发布的 JSON 报告（仅成功任务） |
| `GET` | `/jobs/{job_id}/coverage` | 查询由该质控任务发起的最近一个覆盖任务 |

`min_quality` 为 0..93 的整数（Phred+33 可表示的范围）。上传字段名固定为
`r1`、`r2`。`/coverage` 与 `/coverage/{job_id}` 是 `/coverage-jobs` 的等价
别名，不出现在 Swagger 中。

### 示例

```bash
# 提交任务
curl -X POST http://localhost:8000/jobs \
  -F 'r1=@samples/mixed_R1.fastq.gz' \
  -F 'r2=@samples/mixed_R2.fastq.gz' \
  -F 'min_quality=20'

# 查询
curl http://localhost:8000/jobs/<job_id>

# 停止
curl -X POST http://localhost:8000/jobs/<job_id>/stop

# 成功后可下载
curl -OJ http://localhost:8000/jobs/<job_id>/download/r1
curl -OJ http://localhost:8000/jobs/<job_id>/download/r2

# 基于该成功质控任务提交覆盖核查
curl -X POST http://localhost:8000/coverage-jobs \
  -F 'qc_job_id=<job_id>' \
  -F 'reference=@samples/target_reference.fa' \
  -F 'min_fragment=100' \
  -F 'max_fragment=100' \
  -F 'min_depth=1'

# 查询覆盖任务 / 下载报告
curl http://localhost:8000/coverage-jobs/<coverage_job_id>
curl -OJ http://localhost:8000/coverage-jobs/<coverage_job_id>/report
```

覆盖核查的参考序列也可以作为普通文本表单字段 `reference` 提交（值为 FASTA
文本），gzip 压缩的 FASTA 同样接受。`fragment_length_range=100-300`（或
`100..300`）可替代 `min_fragment`/`max_fragment`。

### 状态字段

质控任务与覆盖任务共用 `status` 状态机：`pending` → `running` →
`succeeded` / `failed`，或 `running`/`pending` → `stopping` → `stopped`。
覆盖任务额外带 `kind: "coverage"`。

成功时：

```json
{
  "status": "succeeded",
  "total_pairs": 7,
  "kept_pairs": 4,
  "discarded_pairs": 3,
  "downloads": {
    "r1": "/jobs/<id>/download/r1",
    "r2": "/jobs/<id>/download/r2"
  },
  "error": null
}
```

失败时 `kept_pairs`/`discarded_pairs`/`downloads` 为 `null`，`error.code`
为下列之一：

| code | 含义 |
|------|------|
| `pair_name_mismatch` | 去掉 `/1`、`/2` 后缀后两端名称不同 |
| `invalid_name_suffix` | 名称未以 `/1`（R1）或 `/2`（R2）结尾 |
| `pair_count_mismatch` | 两端记录条数不同 |
| `seq_qual_length_mismatch` | 序列与质量行长度不同 |
| `malformed_record` | 记录不是完整四行、缺少 `@`/`+`、行超长等 |
| `invalid_quality_char` | 质量字符超出 Phred+33 可打印范围（ASCII 33..126） |
| `gzip_corrupt` | gzip 截断、CRC 错误、尾部垃圾、非 gzip 文件或空文件 |

## 处理规则

1. 每条读段恰好四行：`@名称`、序列、`+`、质量。
2. R1 名称必须以 `/1` 结尾，R2 必须以 `/2` 结尾；去掉后缀后逐对核对，
   必须完全一致。名称后的注释（空格分隔）会被忽略。
3. 序列与质量长度必须一致。
4. 质量按 Phred+33 解释；只要任一端出现低于阈值的碱基，整对丢弃。
5. 保留的读段按原顺序写出，输出名称沿用原始四行（仅补一个 `\n`）。
6. gzip 损坏、条数不符、名称不符、长度不符均导致任务 `failed`，不产生任何
   可下载文件。
7. 流式处理：逐对读取、逐对写出，不会整份载入内存；每 256 对更新一次进度。
   处理期间状态查询与停止请求由独立线程和事件立即响应。
8. 停止与完成互斥：一旦接受停止，任务不可能转为成功；完成后再次停止返回
   `409`。`stopping` 只在工作线程退出且临时文件（上传文件、`.part` 中间
   输出）清空后，才会变为 `stopped`。
9. 只有 `succeeded` 任务提供下载；失败或停止均返回 `409`，不会给出半成品。
   输出先在 `.part` 文件中完整写入，成功后才原子改名发布。

覆盖任务成功时的报告字段（`/coverage-jobs/{id}`，报告文件为其精简版）：

```json
{
  "job_id": "<coverage_id>",
  "kind": "coverage",
  "status": "succeeded",
  "qc_job_id": "<source_qc_id>",
  "reference": {"name": "target_amplicon", "length": 240},
  "min_fragment": 100,
  "max_fragment": 100,
  "min_depth": 1,
  "total_pairs": 3,
  "placed_pairs": 3,
  "ambiguous_pairs": 0,
  "unplaced_pairs": 0,
  "pair_counts": {"placed": 3, "ambiguous": 0, "unplaced": 0},
  "low_coverage_intervals": [{"start": 120, "end": 140}],
  "report": "/coverage-jobs/<id>/report",
  "error": null
}
```

失败/停止的覆盖任务中，上述统计与区间字段均为 `null`。

## 覆盖核查规则

覆盖核查任务只接受**已成功**的质控任务（`/jobs/{id}/download/*` 可用的那些
结果），直接流式重读已发布的 `R1.kept.fastq.gz` / `R2.kept.fastq.gz`，不复制、
不改写质控结果。

**参考序列**：FASTA 必须恰好一条记录，长度 1..50000 bp，仅允许 `ACGT`
（大小写不敏感，行内空白忽略）。多记录、空序列、N 等其它碱基、超过 50 kb
或非 ASCII 内容均返回 `422`。

**落点判定**：逐对处理，R1 按正向、R2 反向互补后在参考上做精确匹配
（不允许错配、不允许多于/少于参考上的完全匹配）。一个“合规组合”需同时满足：

* R1 起点 `p` 不晚于 R2 起点 `q`（`p <= q`）；
* 片段 `[p, q + len(R2))` 完整包含两端读段；
* 片段长度 `q + len(R2) - p` 落在提交的闭区间
  `[min_fragment, max_fragment]` 内。

一对读段：

* 合规组合数 **0** → 记 `unplaced_pairs`（未定位）；
* 合规组合数 **1** → 记 `placed_pairs`，计入覆盖；
* 合规组合数 **≥2** → 记 `ambiguous_pairs`（歧义），**不计入覆盖**。

`placed_pairs + ambiguous_pairs + unplaced_pairs = total_pairs`。

**深度与区间**：只有唯一定位的读段对增加深度。每对把片段
`[R1 起点, R2 末端)` 整段各计一次深度：两端重叠时重叠碱基不会重复计数，
两端之间若有空隙也属于该片段并同样计入。低覆盖区间是参考上深度
**严格小于** `min_depth` 的极大连续区间，坐标 0-based、半开（`start` 含，
`end` 不含）。`min_depth=0` 时没有低覆盖区间。

**失败与停止**：源 gzip 后段损坏（截断、CRC、尾部垃圾）、记录不合法、参考
校验失败或收到停止请求都会使任务 `failed`/`stopped`，**不发布任何报告**；
查询与报告下载返回 `409`，任务目录先清空，之后才报告 `stopped`。报告先写入
`coverage.json.part`，完整成功后才原子改名发布。覆盖核查运行期间，其他质控
任务的查询、下载与停止不受影响。

**失败码**：除质控部分已有的 `gzip_corrupt`、`malformed_record`、
`pair_count_mismatch`、`pair_name_mismatch`、`seq_qual_length_mismatch`
外，覆盖核查还可能出现 `invalid_fragment_range`、`invalid_min_depth`、
`invalid_reference` 与 `internal_error`。

## 本地样例

`tools/make_samples.py` 可重新生成 `samples/` 下的样例：

| 文件 | 内容 |
|------|------|
| `good_R1/R2.fastq.gz` | 3 对全部通过的读段 |
| `bad_R1/R2.fastq.gz` | 3 对全部因低质量被丢弃 |
| `mixed_R1/R2.fastq.gz` | 7 对，前 3 与最后 1 对保留，中间 3 对丢弃（阈值 20） |
| `truncated_R1.fastq.gz` + `truncated_R2.fastq.gz` | R2 的 gzip 流末尾被截断，任务必定失败为 `gzip_corrupt` |
| `target_R1/R2.fastq.gz` + `target_reference.fa` | 240 bp 参考上的 3 对 100 bp 片段；覆盖深度 1 的区域恰好是 `[120, 140)` |

```bash
python tools/make_samples.py
```

对 `target` 样例做覆盖核查（`min_fragment=max_fragment=100`、`min_depth=1`）
应得到：

```json
{
  "status": "succeeded",
  "total_pairs": 3,
  "placed_pairs": 3,
  "ambiguous_pairs": 0,
  "unplaced_pairs": 0,
  "low_coverage_intervals": [{"start": 120, "end": 140}]
}
```

## 测试

```bash
pip install -r requirements.txt pytest httpx
python -m pytest tests -q
```

测试覆盖过滤/顺序/阈值边界、各类强制失败、停止语义（含停止抢先于完成/
失败、清理时序）、处理期间查询响应性、流式读取行数上限，以及覆盖核查的
落点规则（正向/反向互补、重叠只计一次、歧义、未定位、区间过滤）、参考序列
校验、后段 gzip 损坏不发布报告、停止清理时序和查询响应性。

## 说明与限制

* 单进程内存任务表，符合“单个 API 进程管理任务”的范围；不要用多 worker
  启动（`--workers 1`），否则任务查询会落到不同进程。Compose 默认单进程。
* 每个上传文件默认上限 20 GiB，可用 `QC_MAX_UPLOAD_BYTES` 调整；参考序列
  上传上限 8 MiB（解压后同样受限）。
* 后台线程为守护线程；进程关闭时会给所有任务发停止信号并等待清理。
* 覆盖核查直接复用质控任务发布的两份 gzip 结果；清理质控任务目录会使对应
  覆盖任务不可再运行。
