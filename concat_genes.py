#!/usr/bin/env python3
"""
将多个基因的未比对序列按物种串联（每个基因文件内的序列长度可不等，但通常应等长）。
默认要求所有基因文件的物种ID完全一致，否则报错并停止。
用法：
    python concat_unali_strict.py -i 16s.fas rpob.fas ... -o merged_unaligned.fasta
    python concat_unali_strict.py -i *.fas -o merged_unaligned.fasta --allow-missing --missing N
"""

import sys
import argparse
from collections import OrderedDict

def read_fasta(filepath):
    """读取 FASTA 文件，返回 {species: seq}"""
    seqs = OrderedDict()
    current_id = None
    current_seq = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq).upper()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_id is not None:
        seqs[current_id] = ''.join(current_seq).upper()
    return seqs

def main():
    parser = argparse.ArgumentParser(description="串联多个基因的序列（未比对），强制检查物种一致性")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="基因文件列表")
    parser.add_argument("-o", "--output", required=True, help="输出文件")
    parser.add_argument("--allow-missing", action="store_true",
                        help="允许物种缺失，缺失位置用 --missing 字符填充")
    parser.add_argument("--missing", default="", 
                        help="缺失物种的填充字符（默认空字符串，即不填充）")
    args = parser.parse_args()

    input_files = args.input
    print("\n🔗 基因串联顺序：")
    for i, f in enumerate(input_files, 1):
        print(f"   {i}. {f}")

    # 读取所有基因文件
    gene_seqs = []          # 列表元素：{species: seq}
    all_species_sets = []   # 每个文件的物种集合
    gene_lengths = []       # 记录每个基因文件中最长序列的长度（用于缺失填充提示）
    
    for f in input_files:
        seqs = read_fasta(f)
        gene_seqs.append(seqs)
        all_species_sets.append(set(seqs.keys()))
        max_len = max((len(s) for s in seqs.values()), default=0)
        gene_lengths.append(max_len)
        print(f"       ↳ 物种数: {len(seqs)}, 最长序列长度: {max_len} bp")

    # ------------------- 物种 ID 一致性检查 -------------------
    if not args.allow_missing:
        reference_species = all_species_sets[0]
        consistent = True
        for i, sp_set in enumerate(all_species_sets):
            if sp_set != reference_species:
                consistent = False
                missing_in_this = reference_species - sp_set
                extra_in_this = sp_set - reference_species
                if missing_in_this:
                    print(f"\n❌ 错误：文件 '{input_files[i]}' 缺少以下物种（相对于第一个文件）：")
                    for sp in sorted(missing_in_this):
                        print(f"       - {sp}")
                if extra_in_this:
                    print(f"\n❌ 错误：文件 '{input_files[i]}' 多出以下物种（相对于第一个文件）：")
                    for sp in sorted(extra_in_this):
                        print(f"       + {sp}")
        if not consistent:
            sys.exit("\n💥 物种 ID 集合不一致，脚本已终止。如需允许缺失，请添加 --allow-missing 参数。")
    # ---------------------------------------------------------

    # 所有物种的并集（用于允许缺失的情况）
    all_species = set()
    for sp_set in all_species_sets:
        all_species.update(sp_set)
    species_list = sorted(all_species)

    # 串联序列
    concatenated = OrderedDict()
    for sp in species_list:
        parts = []
        for i, seqs in enumerate(gene_seqs):
            if sp in seqs:
                parts.append(seqs[sp])
            else:
                # 缺失物种：填充指定字符或留空
                fill_len = gene_lengths[i] if args.missing else 0
                parts.append(args.missing * fill_len)
        concatenated[sp] = ''.join(parts)

    # 写入文件
    with open(args.output, 'w') as f:
        for sp in species_list:
            f.write(f">{sp}\n")
            seq = concatenated[sp]
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    print(f"\n✅ 串联完成，输出文件: {args.output}")
    print(f"   物种总数: {len(species_list)}")
    print(f"   各基因拼接情况:")
    start = 1
    for fname, max_len in zip(input_files, gene_lengths):
        end = start + max_len - 1 if max_len else start
        print(f"      {fname}: 约 {start}-{end} (基于最长序列)")
        start = end + 1 if max_len else start

    print("\n⚠️  重要提示：")
    print("   此脚本仅简单拼接序列，未检查同源位点对齐。")
    print("   请确保每个基因文件内序列长度相同且为同源片段。")
    print("   建议下一步使用 MAFFT 等工具对输出文件进行多序列比对。")

if __name__ == "__main__":
    main()