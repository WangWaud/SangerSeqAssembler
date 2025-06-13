# Sanger 双向测序结果自动拼接工具

本项目提供一个基于 Python 的命令行工具，可批量自动拼接细菌16S等双向 Sanger 测序（ab1 或 fasta 文件）结果，极大提升实验数据处理效率。

## 功能简介

- 支持 ab1 原始测序文件或 fasta 文件的批量拼接
- 自动识别正向/反向文件并配对
- 自动寻找重叠区并拼接，输出高质量的拼接序列
- 输出拼接质量指标（Overlap、Score）
- 详细命令行帮助文档

## 安装依赖

建议使用 Python 3.7 及以上版本。

只需安装 biopython：

```bash
pip install biopython>=1.81
```

## 使用方法

1. **准备数据文件夹**

   - 将正向测序文件（如 `*_F.ab1` 或 `*_F.fasta`）放入一个文件夹（如 `forward/`）
   - 将反向测序文件（如 `*_R.ab1` 或 `*_R.fasta`）放入另一个文件夹（如 `reverse/`）

2. **运行脚本**

   ```bash
   python sequence_assembly.py --forward_dir ./forward --reverse_dir ./reverse --output_dir ./assembled --file_type .ab1
   ```

   - `--forward_dir`：正向测序文件目录
   - `--reverse_dir`：反向测序文件目录
   - `--output_dir`：拼接结果输出目录
   - `--file_type`：输入文件类型，支持 `.ab1` 或 `.fasta`，默认 `.fasta`

3. **查看帮助文档**

   ```bash
   python sequence_assembly.py -h
   ```

## 文件命名要求

- 正向文件：`样本名_F.ab1` 或 `样本名_F.fasta`
- 反向文件：`样本名_R.ab1` 或 `样本名_R.fasta`
- "样本名"部分需正反向一致，脚本会自动配对
- **示例：**
  - `1_F.ab1` 和 `1_R.ab1`（样本名为 `1`）
  - `A12_F.ab1` 和 `A12_R.ab1`（样本名为 `A12`）
  - `sampleX_F.fasta` 和 `sampleX_R.fasta`（样本名为 `sampleX`）

## 输出说明

- 拼接结果以 `样本名_assembled.fasta` 命名，保存在输出目录
- 每条序列的描述信息包含：
  - **Overlap**：正向序列和反向互补序列之间的重叠碱基数，重叠越长拼接越可靠
  - **Score**：重叠区域的比对得分，反映重叠区的相似度，得分越高说明拼接越可靠

## 典型输出示例

```text
>1_assembled Overlap: 120, Score: 115.0
ATGCCG...
```

## 适用场景

- 细菌16S rDNA全长测序拼接
- 其他需要双向 Sanger 测序拼接的生物信息学场景

## 批量重命名脚本（可选）

如果你的测序文件名为测序公司常见格式（如生工： `编号_..._(样本名)_[引物名].ab1` 或 `.fasta`），可使用本项目自带的 `rename_ab1_files.py` 脚本批量重命名为自动拼接所需格式。

**用法示例：**

```bash
python rename_ab1_files.py --dir ./seqfiles
```

- `--dir`：需要重命名的 ab1 或 fasta 文件所在目录，默认当前目录

重命名后，正向文件为 `样本名_F.ab1` 或 `样本名_F.fasta`，反向文件为 `样本名_R.ab1` 或 `样本名_R.fasta`，便于后续自动拼接。

**命名示例：**
- 输入：`0001_31425061200565_(A12)_[27F].ab1` 和 `0001_31425061200565_(A12)_[1492R].fasta`
- 输出：`A12_F.ab1` 和 `A12_R.fasta`

**命令行帮助文档：**

```text
批量重命名Sanger测序公司返回的ab1或fasta文件，便于自动拼接。

本脚本会将形如：
  编号_..._(样本名)_[引物名].ab1 或 .fasta
自动重命名为：
  样本名_F.ab1 / 样本名_F.fasta（正向，27F）
  样本名_R.ab1 / 样本名_R.fasta（反向，1492R）

用法示例：
  python rename_ab1_files.py --dir ./seqfiles

参数说明：
  --dir    需要重命名的ab1或fasta文件所在目录，默认当前目录

命名示例：
  输入：0001_31425061200565_(A12)_[27F].ab1 和 0001_31425061200565_(A12)_[1492R].fasta
  输出：A12_F.ab1 和 A12_R.fasta
```

## 贡献与反馈

欢迎提交 issue 或 pull request 进行改进！
