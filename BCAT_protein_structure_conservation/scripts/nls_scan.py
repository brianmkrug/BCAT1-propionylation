#!/usr/bin/env python3
"""Rule-based NLS scanning of BCAT1/BCAT2 with positive and negative controls.

Implements the classical, fully published pattern rules:
  * PSORT pat4 / pat7 / bipartite (Nakai & Kanehisa 1992; Nakai & Horton 1999)
  * canonical monopartite consensus K(K/R)X(K/R)  (Chelsky 1989 / Lange 2007)
  * canonical bipartite consensus (K/R)(K/R)-X(9-12)-(K/R)(3 of 5)
    (Dingwall & Laskey 1991; Robbins 1991)
  * sliding-window basic-residue content and net charge
No attempt is made to re-implement cNLS Mapper, NLStradamus or NucPred - those
are run on their own servers and their outputs recorded separately.
"""
import re, json, urllib.request
import numpy as np

BASIC = set("KR")
POS = set("KR"); NEG = set("DE")

TARGETS = {
    "BCAT1":        ("P54687", "target"),
    "BCAT2":        ("O15382", "target"),
    "SV40_largeT":  ("P03070", "positive"),   # PKKKRKV, the textbook monopartite NLS
    "Nucleoplasmin":("P05221", "positive"),   # Xenopus, textbook bipartite NLS
    "c-Myc":        ("P01106", "positive"),   # PAAKRVKLD
    "NPM1_human":   ("P06748", "positive"),
    "CA2_human":    ("P00918", "negative"),   # cytosolic carbonic anhydrase II
    "ALDOA_human":  ("P04075", "negative"),   # cytosolic aldolase A
}


def fetch(acc):
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta",
                                timeout=60) as r:
        t = r.read().decode()
    return "".join(l.strip() for l in t.splitlines() if not l.startswith(">"))


# ------------------------------------------------------------------ PSORT rules
def pat4(s):
    """4 consecutive basics, or 3 basics within a 4-mer that also has H or P."""
    hits = []
    for i in range(len(s) - 3):
        w = s[i:i + 4]
        nb = sum(c in BASIC for c in w)
        if nb == 4:
            hits.append((i + 1, w, "4 basic"))
        elif nb == 3 and any(c in "HP" for c in w):
            hits.append((i + 1, w, "3 basic + H/P"))
    return hits


def pat7(s):
    """P, then within 3 residues a 4-mer containing >=3 basics (PSORT pat7)."""
    hits = []
    for i in range(len(s) - 6):
        if s[i] != "P":
            continue
        w = s[i:i + 7]
        if sum(c in BASIC for c in w[3:7]) >= 3:
            hits.append((i + 1, w, "pat7"))
    return hits


def bipartite(s, spacer=(9, 12), tail=5, need=3):
    """(K/R)(K/R) - X(spacer) - window of `tail` with >= `need` basics."""
    hits = []
    for i in range(len(s) - (2 + spacer[1] + tail)):
        if not (s[i] in BASIC and s[i + 1] in BASIC):
            continue
        for g in range(spacer[0], spacer[1] + 1):
            j = i + 2 + g
            w = s[j:j + tail]
            if len(w) == tail and sum(c in BASIC for c in w) >= need:
                hits.append((i + 1, s[i:j + tail], f"bipartite spacer={g}"))
                break
    return hits


def monopartite_consensus(s):
    """K(K/R)X(K/R) - the classical importin-alpha major-site consensus."""
    return [(m.start() + 1, m.group(), "K(K/R)X(K/R)")
            for m in re.finditer(r"(?=(K[KR].[KR]))", s)
            for _ in [0]] or []


def window_stats(s, w=20):
    n = len(s)
    frac, chg = np.zeros(n), np.zeros(n)
    h = w // 2
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        seg = s[a:b]
        frac[i] = sum(c in BASIC for c in seg) / len(seg)
        chg[i] = (sum(c in POS for c in seg) - sum(c in NEG for c in seg))
    return frac, chg


seqs, report = {}, {}
for name, (acc, role) in TARGETS.items():
    s = fetch(acc)
    seqs[name] = s
    p4, p7, bp = pat4(s), pat7(s), bipartite(s)
    mono = [(m.start() + 1, s[m.start():m.start() + 4]) for m in re.finditer(r"(?=K[KR].[KR])", s)]
    f20, c20 = window_stats(s, 20)
    f10, _ = window_stats(s, 10)
    report[name] = dict(
        acc=acc, role=role, n=len(s),
        pct_basic=100 * sum(c in BASIC for c in s) / len(s),
        pat4=p4, pat7=p7, bipartite=bp, monopartite=mono,
        max_basic_20=float(f20.max()), argmax_20=int(f20.argmax()) + 1,
        max_basic_10=float(f10.max()),
        max_netcharge_20=float(c20.max()),
    )
    print(f"{name:14s} {acc} {role:8s} n={len(s):4d} basic={report[name]['pct_basic']:5.1f}%  "
          f"pat4={len(p4):2d} pat7={len(p7):2d} bipartite={len(bp):2d} "
          f"K(K/R)X(K/R)={len(mono):2d}  max %K+R in 20-aa window={100*f20.max():5.1f}% "
          f"(at {report[name]['argmax_20']})  max net charge(20)={c20.max():+.0f}")

print()
for g in ["BCAT1", "BCAT2"]:
    s = seqs[g]
    print(f"--- {g} N-terminal 1-40: {s[:40]}")
    print(f"    basics in 1-25: {[(i+1, c) for i, c in enumerate(s[:25]) if c in BASIC]}")
    print(f"    acidics in 1-25: {[(i+1, c) for i, c in enumerate(s[:25]) if c in NEG]}")
    net = sum(c in POS for c in s[:25]) - sum(c in NEG for c in s[:25])
    print(f"    net charge 1-25 = {net:+d}; %K+R = {100*sum(c in BASIC for c in s[:25])/25:.1f}%")
    print(f"    all pat4 hits: {report[g]['pat4']}")
    print(f"    all pat7 hits: {report[g]['pat7']}")
    print(f"    all bipartite hits: {report[g]['bipartite']}")
    print(f"    K(K/R)X(K/R) hits: {report[g]['monopartite']}")

json.dump({k: {kk: vv for kk, vv in v.items()} for k, v in report.items()},
          open("/root/bc/data/nls_report.json", "w"), indent=1)
json.dump(seqs, open("/root/bc/data/nls_seqs.json", "w"))
print("\nwrote nls_report.json")
