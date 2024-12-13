#!/usr/bin/env python
"""
Author:        Zhenyu Chen
Created Date:  2022/09/24
Last Modified: 2022/09/25
"""


def get_res_short(res_long):
    res_short = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLN": "Q",
        "GLU": "E",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "HER": "H",
    }
    if res_long in list(res_short.keys()):
        re = res_short[res_long]
    else:
        re = "X"
    return re


def get_seq_(file_path):
    with open(file_path) as f:
        lines = f.read().splitlines()
    header_lines = []
    for line in lines:
        if line.startswith("SEQRES"):
            header_lines.append(line)
    seqs_lines = {}
    chain_id = []
    for header_line in header_lines:
        if header_line[11] not in chain_id:
            chain_id.append(header_line[11])
            seqs_lines[header_line[11]] = []
            seqs_lines[header_line[11]].append(header_line)
        else:
            seqs_lines[header_line[11]].append(header_line)

    seqs = {}
    for key in seqs_lines:
        seqs[key] = ""
        for line in seqs_lines[key]:
            resl = line[18:].split(" ")
            for res in resl:
                res = res.strip()
                if len(res) == 3:
                    seqs[key] += get_res_short(res)
    # seqs_str = ''
    # for key in seqs:
    #     seqs_str+=seqs[key]
    #     seqs_str+=':'
    # seqs_str = seqs_str[:-2]

    return seqs


if __name__ == "__main__":
    """
    Get seqs in pdb files
    Input: pdb_path
    Output:{chain_id:chain_seqs}
    """
    seqs = get_seq_("./header1.pdb")
    print(seqs)
