# hickit：Hi-C 工具集（Rust）

一个用 Rust 编写的快速 Hi-C 工具集：分辨率估计（兼容 Juicer）、.hic（straw）工具，以及对 merged_nodups 的区间过滤。

## 概述

工具集包括：
1. 分辨率（resolution）：从 merged_nodups 或 pairtools .pairs 估计图谱分辨率
2. straw：对 .hic 文件进行 list/dump/effres（observed/NONE/BP）
3. filter：按给定区间提取 merged_nodups 中相关的行

## 特点
- **快速**：比原始的 Bash/awk 实现快 10 - 100 倍
- **内存高效**：使用原子计数器和并行处理
- **灵活**：支持压缩（.gz）和未压缩的输入文件
- **符合标准**：与 Juicer 的 merged_nodups 格式兼容

## 安装

从release中下载预编译版本，或者手动编译，添加到环境变量中

```bash
cargo build --release
```


## 使用方法
使用 merged_nodups 文件的基本用法：
```bash
hickit resolution merged_nodups.txt
```

使用压缩输入：
```bash
hickit resolution merged_nodups.txt.gz
```

从标准输入读取：
```bash
zcat merged_nodups.txt.gz | hickit resolution
```

### 分辨率子命令选项
- `--genome-size <SIZE>`：以碱基对（bp）为单位的总基因组大小（默认值：对于 hg19 为 2428425688）
- `--bin-width <WIDTH>`：以碱基对（bp）为单位的基础 bin 宽度（默认值：50）
- `--prop <PROPORTION>`：良好 bin 的所需比例（默认值：0.8）
- `--count-threshold <COUNT>`：每个 bin 的最小接触次数（默认值：1000）
- `--step-size <SIZE>`：粗略搜索的步长（默认值：1000）
- `--threads <NUM>`：线程数（默认值：自动）

### 示例
```bash
# 人类 hg38 基因组
hickit resolution --genome-size 3137161264 merged_nodups.txt

# 自定义参数
hickit resolution --prop 0.75 --count-threshold 500 merged_nodups.txt

# 使用 8 个线程
hickit resolution --threads 8 merged_nodups.txt.gz
```

### Pairtools .pairs 使用方法
该工具还接受带有标题行（例如 `#chromsize:`）的 pairtools `.pairs` 或 `.pairs.gz` 文件：
```bash
# 直接从.pairs 文件读取
hickit resolution data/mapped.pairs

# 读取压缩的.pairs.gz
hickit resolution data/mapped.pairs.gz
```
- 染色体大小自动从 `.pairs` 标题中推导得出；不需要 `--chrom-size`。
- 作为映射质量的替代指标，仅计算 `pair_type == UU` 的行。
- 注意：自动检测依赖于读取文件路径。如果对 `.pairs` 使用标准输入管道，将跳过标题检测；建议直接传递文件路径。

### 区间过滤 merged_nodups

从 `merged_nodups(.gz)` 中提取与给定区间相关的行（任一端命中即输出到 stdout）。

```bash
# 单段式：CHR:START-END
hickit filter data/merged_nodups.txt ptg000001l:23805-33805 > subset.txt

# 两段式：CHR START-END
hickit filter data/merged_nodups.txt ptg000001l 23805-33805 > subset.txt

# 直接读取 gzip（根据 .gz 扩展名自动解压）
hickit filter data/merged_nodups.txt.gz ptg000001l:23805-33805 > subset.txt

# 从标准输入读取（若为 gzip 请自行解压）
zcat data/merged_nodups.txt.gz | hickit filter - ptg000001l:23805-33805 > subset.txt
```

- `--uniq`：套用与主解析一致的唯一性过滤（要求 `mapq1>0 && mapq2>0` 且 `frag1!=frag2`）。
- 区间为闭区间 `[start, end]`；支持分隔符 `-`、`..` 或 `_`，数字可带千分位逗号（如 `23,805-33,805`）。
- 输出为匹配的原始行，便于下游直接使用。

## 输入格式
该工具期望的 Juicer merged_nodups 格式为制表符分隔的字段：
```
str1  chr1  pos1  frag1  str2  chr2  pos2  frag2  mapq1  cigar1  seq1  mapq2
```

当 `mapq1 > 0`、`mapq2 > 0` 且 `frag1 != frag2` 时计算配对。包括染色体内和染色体间的配对。

也支持带有标题行（例如 `#chromsize:`）的 pairtools `.pairs[.gz]` 格式。对于 `.pairs` 输入，该工具会自动检测标题，从中构建染色体长度，并使用 `chrom1 pos1 chrom2 pos2` 列解析数据行。作为映射质量的替代指标，仅使用 `pair_type == UU` 的行。

## 性能
在具有 8 核和 32GB 内存的典型工作站上：

| 输入大小 | 处理时间 |
|------------|----------------|
| 1 亿个配对 | 约 30 秒 |
| 10 亿个配对 | 约 5 分钟 |

内存使用量与基因组大小呈线性关系（对于人类基因组约为 240MB）。

## 与原始版本的比较
此实现产生的结果与原始 Bash 脚本相同，但在性能上有显著提升：
- **解析**：通过优化的字段提取速度快 5 - 10 倍
- **覆盖构建**：通过并行原子操作速度快 3 - 5 倍
- **分辨率搜索**：通过高效的 bin 聚合速度快 2 - 3 倍

## 许可证
MIT 许可证（与原始 Juicer 实现相同） 
