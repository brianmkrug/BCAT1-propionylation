# SUDIPGXIII (SU-DIPG-XIII) K27M vs KO — motif enrichment (FDR-controlled) & motif-bearing TSSs

Contrast: `S.K27M-S.KO` (RECLU combined-peaks). **Control = K27M, Case = KO**, so
`logFC` is KO/K27M. Motifs labeled **higher_in_K27M** come from the down-in-KO
(= K27M-active) promoter set; **higher_in_KO** from the up-in-KO set. Motifs were
scanned −200/+50 bp around CTSSs (foreground = 200 sequences, background = 3,000).

## 1. Differentially enriched motifs with BH-FDR — `SUDIPGXIII_motifs_FDR.csv`
54 enriched motifs total (46 higher in K27M, 8 higher in KO). Benjamini–Hochberg was
applied across the **full set of 94 discovered/tested motifs** for this contrast
(69 from the K27M-high foreground + 25 from the KO-high foreground). The pipeline's
own binomial enrichment p-values were used as input; a from-scratch binomial
recomputation is included as a cross-check (same significance conclusions — it runs
slightly more extreme than the pipeline value).

All 54 reported motifs pass **q_BH < 0.05**. Strongest: AMD_006 (q=8.9e-21), DREME_027
(q=8.6e-11), DREME_022 (q=1.4e-10), AMD_007 (E-box/ETS-like, q=1.1e-8).

The two TFs of interest both survive FDR:
- **SREBF2** — DREME_031 `ATGGCGTGA`, fg 20/200 vs bg 85/3000, p=1.17e-5, **q_BH=4.2e-5**
- **MYC / MAX::MYC** — DREME_046 `AACACRTGG` (contains the E-box CACGTG), fg 18/200 vs bg 95/3000, p=2.97e-4, **q_BH=7.0e-4**
  (The E-box also appears in AMD_007, which matches the full MYC/MAX/USF/MITF bHLH-Zip family.)

## 2. Top K27M-enriched TSSs carrying the SREBF2 motif — `SUDIPGXIII_K27M_TSS_SREBF2.csv`
All 20 foreground TSSs with the DREME_031 SREBF2 motif; every one is K27M-enriched
(FDR<0.05, higher in K27M), ranked by K27M/KO fold change. Top hits include
**PHYHD1, SLC38A5** (amino-acid transporter), **HAVCR1**, and several unannotated
intergenic CTSSs with very large K27M bias (log2FC ≈ +8–10).

## 3. Top K27M-enriched TSSs carrying the MYC motif — `SUDIPGXIII_K27M_TSS_MYC.csv`
All 18 foreground TSSs with the DREME_046 MYC E-box; top hit is **PAGE5**
(cancer-testis antigen, log2FC(K27M/KO)=+9.9, FDR 4e-96, logCPM 7.9), then
**ARLNC1, MIR4527HG, HBE1**, and multiple unannotated CTSSs.

## Caveats
- FIMO foreground = the K27M-high (down-in-KO) cluster set by construction, so "carries
  motif" TSSs are inherently K27M-biased; ranking is by magnitude of that bias.
- Motif enrichment p-values are single-motif binomial tests from the vendor pipeline;
  BH here controls FDR across the 94-motif family for this one contrast (not across lines).
- TSS coordinates are hg38; gene "." = unannotated cluster.
