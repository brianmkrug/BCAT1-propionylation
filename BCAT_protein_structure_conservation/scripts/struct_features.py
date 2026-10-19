#!/usr/bin/env python3
"""Per-residue structural features for BCAT1/BCAT2 from AlphaFold DB models:
pLDDT, DSSP secondary structure, relative solvent accessibility, and a
sequence-based FoldIndex disorder score."""
import json, subprocess, os

OUT = "/root/bc/data"
ACC = {"BCAT1": "P54687", "BCAT2": "O15382"}

# Kyte-Doolittle hydropathy, scaled 0-1 as in FoldIndex (Prilusky 2005)
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)
CHG = dict(D=-1, E=-1, K=1, R=1)

# DSSP 8-state -> 3-state
SS3 = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E",
       "T": "C", "S": "C", "P": "C", "-": "C", " ": "C", "C": "C"}


def foldindex(seq, win=51):
    """FoldIndex = 2.785<H> - |<R>| - 1.151 ; <0 predicts intrinsic disorder."""
    n = len(seq); h = win // 2; out = []
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        w = seq[a:b]
        H = sum((KD.get(c, 0) + 4.5) / 9.0 for c in w) / len(w)
        R = abs(sum(CHG.get(c, 0) for c in w)) / len(w)
        out.append(2.785 * H - R - 1.151)
    return out


def plddt_from_pdb(path):
    vals, seq3 = {}, {}
    aa3to1 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E",
                  GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F",
                  PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            r = int(line[22:26])
            vals[r] = float(line[60:66])
            seq3[r] = aa3to1.get(line[17:20].strip(), "X")
    n = max(vals)
    return [vals[i] for i in range(1, n + 1)], "".join(seq3[i] for i in range(1, n + 1))


def run_dssp(pdb):
    out = subprocess.run(["mkdssp", "--output-format", "dssp", pdb, "/dev/stdout"],
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[:500])
    ss, acc, started = {}, {}, False
    for line in out.stdout.splitlines():
        if line.startswith("  #  RESIDUE"):
            started = True; continue
        if started and len(line) > 38 and line[13] != "!":
            r = int(line[5:10]); ss[r] = SS3.get(line[16], "C"); acc[r] = float(line[35:38])
    return ss, acc


# Tien et al. 2013 theoretical max ASA
MAXASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224,
              I=197, L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285,
              Y=263, V=174)

res = {}
for g, acc in ACC.items():
    pdb = f"{OUT}/AF_{acc}.pdb"
    plddt, seq = plddt_from_pdb(pdb)
    ref = open(f"{OUT}/{acc}.faa").read().strip()
    assert seq == ref, (g, "model sequence != UniProt")
    ss, sasa = run_dssp(pdb)
    n = len(seq)
    ss_s = "".join(ss.get(i + 1, "C") for i in range(n))
    rsa = [min(1.0, sasa.get(i + 1, 0) / MAXASA.get(seq[i], 200)) for i in range(n)]
    fi = foldindex(seq)
    res[g] = dict(uniprot=acc, seq=seq, plddt=plddt, ss=ss_s, rsa=rsa, foldindex=fi)
    frac = {c: ss_s.count(c) / n for c in "HEC"}
    print(f"{g} ({acc}) n={n}  mean pLDDT={sum(plddt)/n:.1f}  "
          f"H/E/C={frac['H']:.2f}/{frac['E']:.2f}/{frac['C']:.2f}  "
          f"res<70 pLDDT: {sum(1 for v in plddt if v < 70)}  "
          f"FoldIndex<0: {sum(1 for v in fi if v < 0)}")

json.dump(res, open(f"{OUT}/struct_features.json", "w"))
print("wrote struct_features.json")
