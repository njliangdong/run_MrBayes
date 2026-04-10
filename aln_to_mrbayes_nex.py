#!/usr/bin/env python3
"""
Convert FASTA-style .aln (multiple sequence alignment) to NEXUS with MrBayes block.

Example:
  python aln_to_mrbayes_nex.py \
    -i merged.aln -o merged_from_aln.nex \
    --outgroup BACILLUS_VELEZENSIS_BSC16A \
    --ngen 200000 --samplefreq 200 --nchains 4
"""

import argparse
import sys
from collections import OrderedDict
import re

VALID_DNA = set("ACGTURYKMSWBDHVN?-")


def read_fasta(path, allow_duplicates=False):
    seqs = OrderedDict()
    order = []
    warnings = []

    name = None
    buf = []

    def flush(current_name, current_buf, line_no):
        if current_name is None:
            return
        seq = "".join(current_buf).upper()
        if current_name in seqs:
            warnings.append(
                f"Duplicate taxon '{current_name}' near line {line_no}; "
                + ("keeping last" if allow_duplicates else "skipping")
            )
            if allow_duplicates:
                seqs[current_name] = seq
            return
        seqs[current_name] = seq
        order.append(current_name)

    with open(path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush(name, buf, line_no)
                name = line[1:].split()[0].strip()
                buf = []
            else:
                # Remove all whitespace and digits
                line = re.sub(r"[\s\d]", "", line)
                cleaned = "".join(c if c.upper() in VALID_DNA else "N" for c in line)
                buf.append(cleaned)

    flush(name, buf, line_no)

    return seqs, order, warnings


def check_alignment(seqs):
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        print("ERROR: Alignment length mismatch:", file=sys.stderr)
        for k, v in seqs.items():
            print(f"  {k}: {len(v)}", file=sys.stderr)
        sys.exit(1)
    return lengths.pop()


def write_nexus(
    seqs,
    order,
    out_path,
    outgroup=None,
    ngen=100000,
    samplefreq=100,
    nchains=4,
    nst=2,
    rates="invgamma",
    burnin=5000,
):
    ntax = len(seqs)
    nchar = check_alignment(seqs)
    max_name_len = max(len(n) for n in seqs.keys())

    matrix_lines = [f"{name:<{max_name_len}} {seqs[name]}" for name in order]

    content = [
        "#NEXUS",
        "",
        "Begin data;",
        f"    Dimensions ntax={ntax} nchar={nchar};",
        "    Format datatype=dna gap=- missing=?;",
        "    Matrix",
    ]
    content.extend(matrix_lines)
    content.append(";")
    content.append("")
    content.append("BEGIN MRBAYES;")
    if outgroup:
        content.append(f"    outgroup {outgroup};")
    content.append(f"    Lset nst={nst} rates={rates};")
    content.append("    Prset statefreqpr=dirichlet(1,1,1,1);")
    content.append(
        f"    mcmcp savebrlens=yes ngen={ngen} samplefreq={samplefreq} nchains={nchains};"
    )
    content.append("    mcmc;")
    content.append("    sump;")
    content.append(f"    sumt contype=allcompat burnin={burnin};")
    content.append("end;")
    content.append("")

    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(content))


def main():
    parser = argparse.ArgumentParser(
        description="Convert FASTA-style .aln to NEXUS with MrBayes block"
    )
    parser.add_argument("-i", "--input", required=True, help="Input .aln (FASTA) file")
    parser.add_argument("-o", "--output", required=True, help="Output .nex file")
    parser.add_argument("--outgroup", help="Outgroup taxon name")
    parser.add_argument("--ngen", type=int, default=100000, help="MCMC ngen (default 100000)")
    parser.add_argument(
        "--samplefreq", type=int, default=100, help="MCMC samplefreq (default 100)"
    )
    parser.add_argument("--nchains", type=int, default=4, help="MCMC nchains (default 4)")
    parser.add_argument(
        "--nst",
        type=int,
        default=2,
        choices=[1, 2, 6],
        help="Substitution model: 1=JC, 2=HKY, 6=GTR",
    )
    parser.add_argument(
        "--rates",
        default="invgamma",
        choices=["equal", "gamma", "invgamma", "adgamma"],
        help="Rate variation model",
    )
    parser.add_argument(
        "--burnin", type=int, default=5000, help="Burn-in for sumt (default 5000)"
    )
    parser.add_argument(
        "--allow-duplicates",
        action="store_true",
        help="Keep duplicate taxon names (last one wins)",
    )

    args = parser.parse_args()

    seqs, order, warnings = read_fasta(args.input, allow_duplicates=args.allow_duplicates)

    if warnings:
        print("Warnings:", file=sys.stderr)
        for w in warnings:
            print(f"  {w}", file=sys.stderr)

    if len(seqs) < 3:
        sys.exit("ERROR: Need at least 3 taxa.")
    if args.outgroup and args.outgroup not in seqs:
        print(
            f"Warning: outgroup '{args.outgroup}' not found in alignment.",
            file=sys.stderr,
        )

    write_nexus(
        seqs,
        order,
        args.output,
        outgroup=args.outgroup,
        ngen=args.ngen,
        samplefreq=args.samplefreq,
        nchains=args.nchains,
        nst=args.nst,
        rates=args.rates,
        burnin=args.burnin,
    )

    nchar = len(next(iter(seqs.values())))
    print(f"Wrote NEXUS: {args.output} (ntax={len(seqs)}, nchar={nchar})")


if __name__ == "__main__":
    main()
