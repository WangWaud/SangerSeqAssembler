#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
import argparse
import shutil
from collections import Counter

def convert_ab1_to_fasta(ab1_file, output_dir):
    """
    将.ab1文件转换为.fasta文件
    """
    try:
        # 读取.ab1文件
        record = SeqIO.read(ab1_file, "abi")
        
        # 创建输出文件名
        base_name = os.path.splitext(os.path.basename(ab1_file))[0]
        fasta_file = os.path.join(output_dir, f"{base_name}.fasta")
        
        # 保存为fasta格式
        SeqIO.write(record, fasta_file, "fasta")
        return fasta_file
    except Exception as e:
        print(f"转换文件 {ab1_file} 时出错: {str(e)}")
        return None

def find_overlap(seq1, seq2, min_overlap=40, min_score=900, debug=False):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2.0
    aligner.mismatch_score = -5.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -6.0

    best = (0, 0, 0)  # overlap, identity, score
    for overlap in range(min(len(seq1), len(seq2)), min_overlap - 1, -1):
        alignment = aligner.align(seq1[-overlap:], seq2[:overlap])[0]
        seqA = str(seq1[-overlap:])
        seqB = str(seq2[:overlap])
        matches = sum((a == b) and (a != '-') and (b != '-') for a, b in zip(seqA, seqB))
        identity = matches / overlap
        if alignment.score >= min_score and overlap >= min_overlap:
            if alignment.score > best[2]:
                best = (overlap, identity, alignment.score)
        if identity > best[1] or (identity == best[1] and alignment.score > best[2]):
            best = (overlap, identity, alignment.score)
    if best[2] >= min_score and best[0] >= min_overlap:
        if debug:
            print(f"  选用最佳重叠: overlap={best[0]}, identity={best[1]:.2f}, score={best[2]}")
        return best[0], best[2]
    if debug:
        print(f"  未找到满足条件的重叠，最佳尝试：overlap={best[0]}, identity={best[1]:.2f}, score={best[2]}")
    return 0, 0

def assemble_sequences(forward_seq, reverse_seq):
    """
    拼接正向和反向序列，使用Cap3参数标准
    """
    reverse_complement = reverse_seq.reverse_complement()
    overlap, score = find_overlap(forward_seq, reverse_complement)
    if overlap > 0:
        assembled_seq = forward_seq[:-overlap] + reverse_complement
        return assembled_seq, overlap, score
    else:
        return None, 0, 0

def quality_trim(record, min_window=20, min_q=10):
    """
    动态加权质量裁剪：找到首尾连续min_window个碱基Q>=min_q的区段，保留中间高质量区
    """
    seq = str(record.seq)
    quals = record.letter_annotations['phred_quality']
    n = len(seq)
    # 找5'端
    start = 0
    for i in range(n - min_window + 1):
        if min(quals[i:i+min_window]) >= min_q:
            start = i
            break
    # 找3'端
    end = n
    for j in range(n, min_window-1, -1):
        if min(quals[j-min_window:j]) >= min_q:
            end = j
            break
    return seq[start:end], quals[start:end]

def consensus_by_quality(seq1, qual1, seq2, qual2):
    """
    重叠区质量加权共识，仿CAP3论文
    """
    consensus = []
    for b1, q1, b2, q2 in zip(seq1, qual1, seq2, qual2):
        if q1 > q2:
            consensus.append(b1)
        else:
            consensus.append(b2)
    return ''.join(consensus)

def process_file_pair(forward_file, reverse_file, output_dir, sample=None):
    """
    处理一对正向和反向测序文件，支持ab1动态加权裁剪和重叠区质量共识
    """
    # 读取正向
    if forward_file.lower().endswith('.ab1'):
        f_record = SeqIO.read(forward_file, "abi")
        f_seq, f_qual = quality_trim(f_record)
        f_seq = Seq(f_seq)
    else:
        f_seq = SeqIO.read(forward_file, "fasta").seq
        f_qual = [40]*len(f_seq)  # 假定fasta为高质量
    # 读取反向
    if reverse_file.lower().endswith('.ab1'):
        r_record = SeqIO.read(reverse_file, "abi")
        r_seq, r_qual = quality_trim(r_record)
        r_seq = Seq(r_seq)
    else:
        r_seq = SeqIO.read(reverse_file, "fasta").seq
        r_qual = [40]*len(r_seq)
    # 反向互补
    r_seq_rc = r_seq.reverse_complement()
    r_qual_rc = r_qual[::-1]  # 质量值也反转
    # 找重叠
    overlap, score = find_overlap(f_seq, r_seq_rc, debug=False)
    if overlap > 0:
        # 重叠区共识
        f_overlap_seq = f_seq[-overlap:]
        f_overlap_qual = f_qual[-overlap:]
        r_overlap_seq = r_seq_rc[:overlap]
        r_overlap_qual = r_qual_rc[:overlap]
        consensus_overlap = consensus_by_quality(f_overlap_seq, f_overlap_qual, r_overlap_seq, r_overlap_qual)
        assembled_seq = f_seq[:-overlap] + Seq(consensus_overlap) + r_seq_rc[overlap:]
        # 输出
        if sample is not None:
            output_file = os.path.join(output_dir, f"{sample}_assembled.fasta")
        else:
            base_name = os.path.splitext(os.path.basename(forward_file))[0]
            output_file = os.path.join(output_dir, f"{base_name}_assembled.fasta")
        record = SeqRecord(assembled_seq, id=f"{sample}_assembled" if sample else f"{base_name}_assembled", description=f"Overlap: {overlap}, Score: {score}")
        SeqIO.write(record, output_file, "fasta")
        return True, output_file
    else:
        # 只在失败时输出详细调试
        find_overlap(f_seq, r_seq_rc, debug=True)
        print(f"样本{sample}拼接失败，未找到满足条件的重叠区。\n")
        return False, None

def process_directory(input_dir, output_dir, file_type):
    """
    处理目录中的所有文件
    """
    if file_type == '.ab1':
        # 直接返回所有ab1文件的完整路径
        files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.lower().endswith('.ab1')]
        temp_dir = None
        return files, temp_dir
    # 下面是fasta的处理逻辑
    temp_dir = os.path.join(output_dir, "temp_fasta")
    os.makedirs(temp_dir, exist_ok=True)
    files = [f for f in os.listdir(input_dir) if f.endswith(file_type)]
    converted_files = []
    for file in files:
        # 转换.ab1文件为.fasta
        if file_type == '.ab1':
            fasta_file = convert_ab1_to_fasta(os.path.join(input_dir, file), temp_dir)
            if fasta_file:
                converted_files.append(fasta_file)
        else:
            # 复制.fasta文件到临时目录
            src = os.path.join(input_dir, file)
            dst = os.path.join(temp_dir, file)
            shutil.copy2(src, dst)
            converted_files.append(dst)
    return converted_files, temp_dir

def main():
    parser = argparse.ArgumentParser(
        description='''\
双向Sanger测序结果自动拼接工具

本工具可批量处理正向和反向的ab1或fasta测序文件，自动识别样本名并拼接，输出拼接后的fasta文件。

使用示例：
  python sequence_assembly.py --forward_dir ./forward --reverse_dir ./reverse --output_dir ./assembled --file_type .ab1

参数说明：
  --forward_dir   正向测序文件目录（如 *_F.ab1 或 *_F.fasta），默认 ./forward
  --reverse_dir   反向测序文件目录（如 *_R.ab1 或 *_R.fasta），默认 ./reverse
  --output_dir    拼接结果输出目录
  --file_type     输入文件类型，支持 .ab1 或 .fasta，默认 .fasta

文件命名要求：
  正向文件：样本名_F.ab1 或样本名_F.fasta
  反向文件：样本名_R.ab1 或样本名_R.fasta
  （样本名部分需正反向一致）

拼接结果说明：
  输出的fasta文件描述中包含：
    Overlap：正向序列和反向互补序列之间的重叠碱基数，重叠越长拼接越可靠。
    Score  ：重叠区域的比对得分，反映重叠区的相似度，得分越高说明拼接越可靠。
''',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--forward_dir', default='./forward', help='正向测序文件目录（如 *_F.ab1 或 *_F.fasta），默认 ./forward')
    parser.add_argument('--reverse_dir', default='./reverse', help='反向测序文件目录（如 *_R.ab1 或 *_R.fasta），默认 ./reverse')
    parser.add_argument('--output_dir', required=True, help='拼接结果输出目录')
    parser.add_argument('--file_type', choices=['.ab1', '.fasta'], default='.fasta',
                      help='输入文件类型，支持 .ab1 或 .fasta，默认 .fasta')
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 处理正向和反向文件
    forward_files, forward_temp_dir = process_directory(args.forward_dir, args.output_dir, args.file_type)
    reverse_files, reverse_temp_dir = process_directory(args.reverse_dir, args.output_dir, args.file_type)
    
    # 获取所有正向和反向文件的样本名
    forward_dict = {}
    print('正向文件提取到的样本名:')
    for f in forward_files:
        basename = os.path.basename(f)
        if basename.lower().endswith('_f.ab1') or basename.lower().endswith('_f.fasta'):
            sample = basename.rsplit('_F.', 1)[0].strip()
            print(f'  {basename} -> "{sample}"')
            forward_dict[sample] = f
    reverse_dict = {}
    print('反向文件提取到的样本名:')
    for r in reverse_files:
        basename = os.path.basename(r)
        if basename.lower().endswith('_r.ab1') or basename.lower().endswith('_r.fasta'):
            sample = basename.rsplit('_R.', 1)[0].strip()
            print(f'  {basename} -> "{sample}"')
            reverse_dict[sample] = r

    # 找到正反向都存在的样本
    common_samples = set(forward_dict.keys()) & set(reverse_dict.keys())
    print(f'正反向共有样本名: {sorted(common_samples)}')
    if not common_samples:
        print('未找到任何正反向都存在的样本，检查文件命名和目录！')
    for sample in sorted(common_samples):
        forward_file = forward_dict[sample]
        reverse_file = reverse_dict[sample]
        print(f'正在处理样本: {sample} -> {os.path.basename(forward_file)}, {os.path.basename(reverse_file)}')
        success, output_file = process_file_pair(forward_file, reverse_file, args.output_dir, sample)
        if success:
            print(f"成功拼接: {os.path.basename(forward_file)} 和 {os.path.basename(reverse_file)}")
            print(f"输出文件: {output_file}")
        else:
            print(f"拼接失败: {os.path.basename(forward_file)} 和 {os.path.basename(reverse_file)}")
    # 清理临时目录
    if forward_temp_dir is not None and os.path.exists(forward_temp_dir):
        shutil.rmtree(forward_temp_dir)
    if reverse_temp_dir is not None and os.path.exists(reverse_temp_dir):
        shutil.rmtree(reverse_temp_dir)

if __name__ == "__main__":
    main() 
