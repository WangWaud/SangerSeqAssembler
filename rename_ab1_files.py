import os
import re
import argparse

def main():
    parser = argparse.ArgumentParser(
        description='''\
批量重命名Sanger测序公司返回的ab1文件，便于自动拼接。

本脚本会将形如：
  编号_..._(样本名)_[引物名].ab1
自动重命名为：
  样本名_F.ab1（正向，27F）
  样本名_R.ab1（反向，1492R）

用法示例：
  python rename_ab1_files.py --dir ./ab1files

参数说明：
  --dir    需要重命名的ab1文件所在目录，默认当前目录

命名示例：
  输入：0001_31425061200565_(A12)_[27F].ab1 和 0001_31425061200565_(A12)_[1492R].ab1
  输出：A12_F.ab1 和 A12_R.ab1
''',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--dir', default='.', help='ab1文件所在目录，默认当前目录')
    args = parser.parse_args()

    folder = args.dir
    for filename in os.listdir(folder):
        if filename.endswith('.ab1'):
            # 匹配 (样本名) 和引物名
            match = re.match(r'.*\(([^)]+)\)_\[(27F|1492R)\]\.ab1', filename)
            if match:
                sample = match.group(1)
                primer = match.group(2)
                if primer == '27F':
                    newname = f'{sample}_F.ab1'
                elif primer == '1492R':
                    newname = f'{sample}_R.ab1'
                else:
                    continue
                os.rename(os.path.join(folder, filename), os.path.join(folder, newname))
                print(f'{filename} -> {newname}')

if __name__ == '__main__':
    main() 