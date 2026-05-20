#!/usr/bin/env python3
"""
Author : lb21z855 <lb21z855@ifik-srvngs02.ifik.unibe.ch>
Date   : 2024-06-24
Purpose: generate sequences
"""

import argparse
from random import choices


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="generate sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--flank",
        help="Size of the flanking region missing in the query read",
        metavar="int",
        type=int,
        default=20,
    )

    parser.add_argument(
        "-s",
        "--size",
        help="Size of the sequence",
        metavar="int",
        type=int,
        default=150,
    )

    parser.add_argument(
        "-f",
        "--file",
        help="the output file",
        metavar="FILE",
        type=argparse.FileType("wt"),
        default="seq.fasta",
    )

    parser.add_argument(
        "-r",
        "--reverse",
        help="Return reve complement of seq as well",
        action="store_true",
    )

    return parser.parse_args()


def RandSeq(size: int):
    return "".join(choices("ATCG", k=size))


def RevCom(seq: str):
    seq_inv = seq[::-1]
    bases_dict: dict[str, str] = {"A": "T", "C": "G", "T": "A", "G": "C"}
    return "".join([bases_dict[nuc] for nuc in list(seq_inv)])


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    flank = args.flank
    # size = args.size

    with args.file as fasta:
        seq: str = RandSeq(size=args.size)
        rev_seq_1: str = RevCom(seq=seq)  # case1
        rev_seq_2: str = RevCom(seq[flank:])  # case2
        rev_seq_3: str = RevCom(seq[flank:-flank])  # case 3
        rev_seq_4: str = RevCom(seq[:-flank])  # case 4

        seq_1: str = seq
        seq_2: str = seq[flank:]
        seq_3: str = seq[flank:-flank]
        seq_4: str = seq[:-flank]

        if args.reverse:
            print(f">seq_1\n{seq_1}", file=fasta)
            print(f">seq_2\n{seq_2}", file=fasta)
            print(f">seq_3\n{seq_3}", file=fasta)
            print(f">seq_4\n{seq_4}", file=fasta)

            print(f">rev_seq_1\n{rev_seq_1}", file=fasta)
            print(f">rev_seq_2\n{rev_seq_2}", file=fasta)
            print(f">rev_seq_3\n{rev_seq_3}", file=fasta)
            print(f">rev_seq_4\n{rev_seq_4}", file=fasta)
        else:
            print(f">seq_1\n{seq_1}", file=fasta)
            print(f">seq_2\n{seq_2}", file=fasta)
            print(f">seq_3\n{seq_3}", file=fasta)
            print(f">seq_4\n{seq_4}", file=fasta)


# --------------------------------------------------
if __name__ == "__main__":
    main()
