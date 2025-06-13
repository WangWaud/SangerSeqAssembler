import os
import re
import argparse
import shutil

def main():
    parser = argparse.ArgumentParser(
        description='''\
批量重命名并分拣Sanger测序公司返回的ab1或fasta文件，便于自动拼接。

本脚本会将形如：
  编号_..._(样本名)_[引物名].ab1 或 .fasta
自动重命名为：
  样本名_F.ab1 / 样本名_F.fasta（正向，27F）
  样本名_R.ab1 / 样本名_R.fasta（反向，1492R）
并自动将正向文件移动到forward/，反向文件移动到reverse/。

用法示例：
  python rename_ab1_files.py --dir ./seqfiles

参数说明：
  --dir    需要重命名和分拣的ab1或fasta文件所在目录，默认当前目录

命名示例：
  输入：0001_31425061200565_(A12)_[27F].ab1 和 0001_31425061200565_(A12)_[1492R].fasta
  输出：forward/A12_F.ab1 和 reverse/A12_R.fasta
''',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--dir', default='.', help='ab1或fasta文件所在目录，默认当前目录')
    args = parser.parse_args()

    folder = args.dir
    forward_dir = os.path.join(folder, 'forward')
    reverse_dir = os.path.join(folder, 'reverse')
    os.makedirs(forward_dir, exist_ok=True)
    os.makedirs(reverse_dir, exist_ok=True)

    for filename in os.listdir(folder):
        if filename.endswith('.ab1') or filename.endswith('.fasta'):
            filepath = os.path.join(folder, filename)
            if os.path.isdir(filepath):
                continue
            # 匹配 (样本名) 和引物名
            match = re.match(r'.*\(([^)]+)\)_\[(27F|1492R)\]\.(ab1|fasta)', filename)
            if match:
                sample = match.group(1)
                primer = match.group(2)
                ext = match.group(3)
                if primer == '27F':
                    newname = f'{sample}_F.{ext}'
                    target_dir = forward_dir
                elif primer == '1492R':
                    newname = f'{sample}_R.{ext}'
                    target_dir = reverse_dir
                else:
                    continue
                target_path = os.path.join(target_dir, newname)
                shutil.move(filepath, target_path)
                print(f'{filename} -> {os.path.relpath(target_path, folder)}')

if __name__ == '__main__':
    main() 
