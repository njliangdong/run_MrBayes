#!/usr/bin/env python3
"""
清理 FASTA 文件序列 ID 和非法序列字符。
- ID 清理：仅保留字母、数字、下划线，其他字符替换为下划线。保留 > 之后的全部内容。
- 序列清理：根据指定模式（DNA/RNA/Protein）将非法字符替换为 N（或 X）。
- 默认跳过重复序列（保留第一条），可选用 --rename-duplicates 重命名。
- 支持通配符批量处理（如 *.fasta）。
- 详细报告 ID 修改情况和序列非法字符统计。

用法：
    python clean_fasta_ids.py [--seqtype {dna,rna,protein}] [--rename-duplicates] file1.fasta *.fas ...
"""

import sys
import re
import glob
from collections import OrderedDict

VALID_CHARS = {
    'dna': set("ACGTURYKMSWBDHVN?-"),
    'rna': set("ACGURYKMSWBDHVN?-"),
    'protein': set("ACDEFGHIKLMNPQRSTVWY*-?")
}

REPLACE_CHAR = {
    'dna': 'N',
    'rna': 'N',
    'protein': 'X'
}

def clean_id(raw_id):
    """将非法字符替换为下划线，只保留字母数字下划线"""
    return re.sub(r'[^A-Za-z0-9_]', '_', raw_id)

def clean_sequence(seq_str, seq_type):
    allowed = VALID_CHARS[seq_type]
    replace_char = REPLACE_CHAR[seq_type]
    cleaned = []
    replaced = 0
    for ch in seq_str.upper():
        if ch in allowed:
            cleaned.append(ch)
        else:
            cleaned.append(replace_char)
            replaced += 1
    return ''.join(cleaned), replaced

def process_fasta(input_path, output_path=None, rename_duplicates=False, seq_type='dna'):
    if output_path is None:
        output_path = re.sub(r'\.(fasta|fas|fa)$', '_clean.fas', input_path)
        if output_path == input_path:
            output_path = input_path + '_clean.fas'
    
    seqs = OrderedDict()
    order = []
    warnings = []
    id_changes = []
    current_id = None
    current_raw_id = None
    current_seq_parts = []
    total_original = 0
    total_replaced = 0
    
    with open(input_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                total_original += 1
                # 处理上一条序列
                if current_id is not None:
                    full_seq = ''.join(current_seq_parts)
                    clean_seq, replaced = clean_sequence(full_seq, seq_type)
                    total_replaced += replaced
                    if current_id in seqs:
                        if rename_duplicates:
                            base = current_id
                            count = 1
                            while current_id in seqs:
                                current_id = f"{base}_{count}"
                                count += 1
                            warnings.append(f"名称 '{base}' 重复，已重命名为 '{current_id}'")
                            id_changes.append((current_raw_id, current_id, "重命名", replaced))
                            seqs[current_id] = clean_seq
                            order.append(current_id)
                        else:
                            warnings.append(f"名称 '{current_id}' 重复出现（第 {line_num} 行），已跳过")
                            id_changes.append((current_raw_id, current_id, "跳过（重复）", replaced))
                    else:
                        seqs[current_id] = clean_seq
                        order.append(current_id)
                        action = "修改" if current_raw_id != current_id else "未修改"
                        id_changes.append((current_raw_id, current_id, action, replaced))
                
                # 读取新序列头：保留 > 之后的全部内容（修复点）
                raw_id = line[1:].strip()
                current_raw_id = raw_id
                current_id = clean_id(raw_id)
                current_seq_parts = []
            else:
                current_seq_parts.append(line.upper())
    
    # 处理最后一条序列（同上）
    if current_id is not None:
        full_seq = ''.join(current_seq_parts)
        clean_seq, replaced = clean_sequence(full_seq, seq_type)
        total_replaced += replaced
        if current_id in seqs:
            if rename_duplicates:
                base = current_id
                count = 1
                while current_id in seqs:
                    current_id = f"{base}_{count}"
                    count += 1
                warnings.append(f"名称 '{base}' 重复，已重命名为 '{current_id}'")
                id_changes.append((current_raw_id, current_id, "重命名", replaced))
                seqs[current_id] = clean_seq
                order.append(current_id)
            else:
                warnings.append(f"名称 '{current_id}' 重复出现（文件末尾），已跳过")
                id_changes.append((current_raw_id, current_id, "跳过（重复）", replaced))
        else:
            seqs[current_id] = clean_seq
            order.append(current_id)
            action = "修改" if current_raw_id != current_id else "未修改"
            id_changes.append((current_raw_id, current_id, action, replaced))
    
    # 写入输出文件
    with open(output_path, 'w') as f:
        for sid in order:
            f.write(f">{sid}\n")
            seq = seqs[sid]
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
    
    # 打印详细报告（显示时适当截断长ID）
    print(f"\n📄 文件: {input_path} -> {output_path}")
    print(f"🧬 序列类型: {seq_type.upper()}，非法字符替换为 '{REPLACE_CHAR[seq_type]}'")
    print(f"🔧 ID 校正报告 (共 {len(id_changes)} 条序列):")
    print(f"   {'原始 ID':<40} {'清理后 ID':<40} {'状态':<12} {'替换数'}")
    print(f"   {'-'*40} {'-'*40} {'-'*12} {'-'*6}")
    for raw, clean, action, rep in id_changes:
        raw_short = raw if len(raw) <= 40 else raw[:37] + "..."
        clean_short = clean if len(clean) <= 40 else clean[:37] + "..."
        print(f"   {raw_short:<40} {clean_short:<40} {action:<12} {rep}")
    
    if warnings:
        print(f"⚠️  重复警告:")
        for w in warnings:
            print(f"   {w}")
    
    print(f"✅ 处理完成: 原始序列 {total_original} 条，保留 {len(order)} 条")
    print(f"   总共替换非法序列字符: {total_replaced} 个")

def print_help():
    help_text = """
用法: python clean_fasta_ids.py [选项] <文件或通配符...>

选项:
  --seqtype {dna,rna,protein}   指定序列类型，默认为 dna。
  --rename-duplicates           遇到重复序列 ID 时自动重命名（添加 _1, _2 后缀）。
  -h, --help                    显示本帮助信息。

输入文件:
  支持直接指定文件名或使用通配符，例如:
    python clean_fasta_ids.py 16s.fasta
    python clean_fasta_ids.py *.fas

输出文件:
  默认输出为原文件名去除 .fasta/.fas/.fa 后加 _clean.fas。

功能说明:
  1. 序列 ID 清理：保留 > 之后的全部内容，将非法字符替换为下划线。
  2. 序列字符清理：根据类型将非法碱基/氨基酸替换为 N 或 X。
  3. 重复序列处理：默认跳过重复 ID 并警告，可重命名保留。
"""
    print(help_text)

def main():
    if '-h' in sys.argv or '--help' in sys.argv:
        print_help()
        sys.exit(0)
    
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    rename_mode = False
    seq_type = 'dna'
    patterns = []
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == "--rename-duplicates":
            rename_mode = True
            i += 1
        elif arg == "--seqtype":
            if i+1 >= len(sys.argv):
                print("错误：--seqtype 需要指定类型 (dna/rna/protein)", file=sys.stderr)
                sys.exit(1)
            typ = sys.argv[i+1].lower()
            if typ not in ('dna', 'rna', 'protein'):
                print("错误：--seqtype 必须是 dna, rna 或 protein", file=sys.stderr)
                sys.exit(1)
            seq_type = typ
            i += 2
        else:
            patterns.append(arg)
            i += 1
    
    files = []
    for pat in patterns:
        matched = glob.glob(pat)
        if not matched:
            print(f"⚠️  警告: 未找到匹配 '{pat}' 的文件，已跳过", file=sys.stderr)
        else:
            files.extend(matched)
    
    if not files:
        print("❌ 错误：没有找到任何要处理的 FASTA 文件。", file=sys.stderr)
        sys.exit(1)
    
    seen = set()
    unique_files = []
    for f in files:
        if f not in seen:
            seen.add(f)
            unique_files.append(f)
    
    print(f"📁 将处理 {len(unique_files)} 个文件，模式: {seq_type.upper()}")
    for f in unique_files:
        print(f"   {f}")
    
    for f in unique_files:
        process_fasta(f, rename_duplicates=rename_mode, seq_type=seq_type)

if __name__ == "__main__":
    main()