# Sanger 双向测序结果自动拼接工具

本项目提供一个基于 Python 的命令行工具，可批量自动拼接细菌16S等双向 Sanger 测序（ab1 或 fasta 文件）结果，极大提升实验数据处理效率。

## 功能简介

- 支持 ab1 原始测序文件或 fasta 文件的批量拼接
- 自动识别正向/反向文件并配对
- 自动寻找重叠区并拼接，输出高质量的拼接序列
- 输出拼接质量指标（Overlap、Score）
- ab1文件自动进行首尾低质量区裁剪，拼接更准确
- 拼接算法严格参照CAP3/SnapGene参数，采用全局最优score拼接点

## 安装依赖

建议使用 Python 3.7 及以上版本。

只需安装 biopython：

```bash
pip install biopython>=1.81
```

## 使用方法

**1. 样本重命名**

在使用拼接脚本前，你需要确保正向和反向文件的命名规范一致。你可以选择**手动重命名**或使用本项目自带的**批量重命名脚本**：

### 方式一：手动重命名

- 将正向文件命名为 `样本名_F.ab1` 或 `样本名_F.fasta`
- 将反向文件命名为 `样本名_R.ab1` 或 `样本名_R.fasta`
- "样本名"部分需正反向一致
- 手动将正向文件放入 `forward/`，反向文件放入 `reverse/`

### 方式二：批量重命名并分拣（推荐）

- 直接运行 `rename_seq_files.py`，会自动将正向文件（如 `*_F.ab1` 或 `*_F.fasta`）移动到 `forward/`，反向文件（如 `*_R.ab1` 或 `*_R.fasta`）移动到 `reverse/`，无需手动分拣。

```bash
# python rename_seq_files.py -h 查看帮助文档
python rename_seq_files.py --dir ./seqfiles
```
- `--dir`：需要重命名和分拣的 ab1 或 fasta 文件所在目录，默认当前目录

命名示例：
- 输入：`0001_31425061200565_(A12)_[27F].ab1` 和 `0001_31425061200565_(A12)_[1492R].fasta`
- 输出：`forward/A12_F.ab1` 和 `reverse/A12_R.fasta`

**rename_seq_files.py 帮助文档：**

```python
批量重命名并分拣Sanger测序公司返回的ab1或fasta文件，便于自动拼接。

本脚本会将形如：
  编号_..._(样本名)_[引物名].ab1 或 .fasta
自动重命名为：
  样本名_F.ab1 / 样本名_F.fasta（正向，27F）
  样本名_R.ab1 / 样本名_R.fasta（反向，1492R）
并自动将正向文件移动到forward/，反向文件移动到reverse/。

用法示例：
  python rename_seq_files.py --dir ./seqfiles

参数说明：
  --dir    需要重命名和分拣的ab1或fasta文件所在目录，默认当前目录

命名示例：
  输入：0001_31425061200565_(A12)_[27F].ab1 和 0001_31425061200565_(A12)_[1492R].fasta
  输出：forward/A12_F.ab1 和 reverse/A12_R.fasta
```

---

**2. 运行拼接脚本**

   ```bash
   # python sequence_assembly.py -h 查看帮助文档
   python sequence_assembly.py --output_dir ./assembled --file_type .ab1
   ```

   - `--forward_dir`：正向测序文件目录，默认 `./forward`
   - `--reverse_dir`：反向测序文件目录，默认 `./reverse`
   - `--output_dir`：拼接结果输出目录（必填）
   - `--file_type`：输入文件类型，支持 `.ab1` 或 `.fasta`，推荐使用 `.ab1` 以获得更高质量拼接

## 拼接算法说明（核心原理）

- **ab1文件自动质量裁剪**：自动识别并裁剪5'和3'端低质量区段（如窗口20bp内Q≥10），只保留高质量区段用于拼接。
- **全局最优score拼接点**：遍历所有可能的重叠长度（overlap≥40），对每个overlap用CAP3参数（match=2, mismatch=-5, gap=-6）比对，选取score最高且score≥900的重叠区作为拼接点。
- **重叠区共识生成**：在最佳重叠区内，对每个位点，直接选择质量值更高的那一端的碱基作为拼接结果（不输出N）。
- **拼接结果**：正向非重叠区 + 重叠区共识 + 反向非重叠区，生成的fasta文件命名为`样本名_assembled.fasta`，描述中包含重叠长度和score。
- **拼接失败时**：输出最佳尝试的overlap、identity和score，便于调试。

## 输出说明

- 拼接结果以 `样本名_assembled.fasta` 命名，保存在输出目录
- 每条序列的描述信息包含：
  - **Overlap**：正向序列和反向互补序列之间的重叠碱基数，重叠越长拼接越可靠
  - **Score**：重叠区域的比对得分，反映重叠区的相似度，得分越高说明拼接越可靠

## 典型输出示例

```text
>1_assembled Overlap: 665, Score: 1288.0
CGAGCCCTTCGGGGTTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCTTGGGTT...
```

## 推荐使用建议

- **优先使用 ab1 文件**，可获得自动质量裁剪和更高准确度。
- 若只有 fasta 文件，也可拼接，但无法自动裁剪低质量区。
- 拼接参数严格参照CAP3/SnapGene，结果高度一致。

## 贡献与反馈

欢迎提交 issue 或 pull request 进行改进！
