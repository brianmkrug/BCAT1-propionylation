#!/usr/bin/env python3
"""Pull MANE CDS structure + UCSC 100-way phyloP/phastCons for BCAT1 and BCAT2,
map every coding base to a codon/protein residue, and verify the translation."""
import json, urllib.request, time, os

API = "https://api.genome.ucsc.edu"
OUT = "/root/bc/data"

GENES = {
    "BCAT1": dict(chrom="chr12", start=24810023, end=24949101, uniprot="P54687"),
    "BCAT2": dict(chrom="chr19", start=48795064, end=48811029, uniprot="O15382"),
}

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in _b:
    for b in _b:
        for c in _b:
            CODON[a + b + c] = _aa[i]; i += 1
COMP = str.maketrans("ACGTN", "TGCAN")


def get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))


def mane_record(g, info):
    d = get(f"{API}/getData/track?genome=hg38;track=mane;chrom={info['chrom']};"
            f"start={info['start']};end={info['end']}")
    recs = [r for r in d["mane"] if r.get("geneName2") == g]
    assert len(recs) == 1, (g, len(recs))
    return recs[0]


def cds_blocks(r):
    """Genomic CDS intervals (0-based half-open), ascending."""
    starts = [int(x) for x in r["chromStarts"].rstrip(",").split(",")]
    sizes = [int(x) for x in r["blockSizes"].rstrip(",").split(",")]
    cs, ce = r["thickStart"], r["thickEnd"]
    out = []
    for s, sz in zip(starts, sizes):
        a = r["chromStart"] + s
        b = a + sz
        a, b = max(a, cs), min(b, ce)
        if b > a:
            out.append((a, b))
    return out


def seq(chrom, s, e):
    d = get(f"{API}/getData/sequence?genome=hg38;chrom={chrom};start={s};end={e}")
    return d["dna"].upper()


def track(name, chrom, s, e):
    d = get(f"{API}/getData/track?genome=hg38;track={name};chrom={chrom};start={s};end={e}")
    vals = {}
    for row in d.get(name, []):
        v = row["value"]
        for p in range(row["start"], row["end"]):
            vals[p] = v
    return vals


results = {}
for g, info in GENES.items():
    r = mane_record(g, info)
    blocks = cds_blocks(r)
    strand = r["strand"]
    chrom = r["chrom"]
    cds_len = sum(b - a for a, b in blocks)
    print(f"{g}: {chrom}:{r['chromStart']}-{r['chromEnd']} {strand} "
          f"{r['ncbiId']} exons={r['blockCount']} CDSbases={cds_len}")

    # genomic positions in transcript order
    pos = []
    for a, b in blocks:
        pos.extend(range(a, b))
    if strand == "-":
        pos = pos[::-1]

    # sequence in transcript order
    dna = ""
    for a, b in blocks:
        dna += seq(chrom, a, b)
    if strand == "-":
        dna = dna.translate(COMP)[::-1]
    prot = "".join(CODON.get(dna[i:i+3], "X") for i in range(0, len(dna) - 2, 3))

    up = "".join(l.strip() for l in open(f"{OUT}/{info['uniprot']}.faa"))
    ok = prot.rstrip("*") == up
    print(f"   translated {len(prot)} aa  stop={prot.endswith('*')}  matches UniProt: {ok}")
    if not ok:
        for i, (x, y) in enumerate(zip(prot, up)):
            if x != y:
                print("   first mismatch at", i + 1, x, y); break

    pp, pc = {}, {}
    for a, b in blocks:
        pp.update(track("phyloP100way", chrom, a, b))
        pc.update(track("phastCons100way", chrom, a, b))
    print(f"   phyloP covered {len(pp)}/{cds_len}  phastCons {len(pc)}/{cds_len}")

    results[g] = dict(
        chrom=chrom, strand=strand, blocks=blocks, transcript=r["ncbiId"],
        ensembl=r["name"], protein=prot.rstrip("*"), uniprot_match=ok,
        pos=pos,
        phyloP=[pp.get(p) for p in pos],
        phastCons=[pc.get(p) for p in pos],
    )

json.dump(results, open(f"{OUT}/conservation.json", "w"))
print("wrote conservation.json")
