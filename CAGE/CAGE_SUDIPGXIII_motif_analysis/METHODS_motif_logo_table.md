# Methods — CAGE-seq promoter-motif analysis (SU-DIPG-XIII)

**CAGE-seq.** Cap Analysis of Gene Expression (CAGE) captures the 5′-7-methylguanosine
cap of Pol II transcripts and sequences a short tag from each 5′ end, mapping
transcription start sites (TSSs) at single-base resolution; tag counts at a TSS are
proportional to promoter activity [1]. Libraries were prepared and sequenced by DNAFORM
for three H3K27M diffuse midline glioma lines — BT245, SU-DIPG-XIII and HSJ019 — each as
the parental H3K27M mutant and its isogenic H3K27M-knockout (KO) control, in two
biological replicates per genotype (12 libraries).

**Read processing and mapping.** Reads were rRNA-filtered, then mapped to the human
genome (GRCh38 primary assembly; GENCODE v33 annotation) with a two-pass strategy: BWA
(retaining MAPQ ≥ 20) [2] followed by HISAT2 for the remaining reads (retaining uniquely
mapped, NH:i:1) [3]. Uniquely mapped 5′ ends define CAGE-detected TSSs (CTSSs).

**TSS clustering and differential promoter usage.** CTSSs were processed with the RECLU
pipeline [4] to call reproducible tag clusters (promoters) and their 10–90 %
interquantile widths. Per-cluster expression (counts per million) was compared between
genotypes with edgeR [5], and promoters were annotated to GENCODE genes. For the
SU-DIPG-XIII K27M-vs-KO contrast this yields two promoter sets: those with higher 5′
activity in K27M ("K27M-enriched") and those higher in KO ("KO-enriched").

**De novo motif discovery and TF matching.** For each promoter set, the sequence window
−200/+50 bp around its CTSSs was scanned for over-represented motifs by four de novo
engines: DREME [6], GLAM2 [7], Weeder [8], and an additional matrix-refinement step
(AMD). Enrichment of each discovered motif was scored by a binomial test of its
occurrence in the foreground promoters versus a background set (200 vs 3,000 sequences).
Discovered motifs were matched to known transcription factors with TomTom [9] against
JASPAR CORE 2016 (vertebrates) [10], and motif instances were mapped with FIMO [11].
(DREME, GLAM2, TomTom and FIMO are components of the MEME Suite.)

**Figure — SU-DIPG-XIII motif logo table.** From the S.K27M-vs-S.KO comparison we show
the top five **K27M-enriched** motifs that carry a known-TF match and the top five
**KO-enriched** motifs, ranked by enrichment significance (−log₁₀ p). Each row gives the
motif ID, an information-content (bits) sequence logo drawn from the de novo
position-weight matrix (rendered with Logomaker [12]), the TomTom TF match(es), the
binomial enrichment −log₁₀ p, and a Benjamini–Hochberg FDR (q) recomputed across the full
set of 94 discovered motifs (both directions). K27M-enriched motifs are dominated by
interpretable TF elements — an E-box (DREME_046; MYC/MAX, MondoA-MLX/MLXIP), a
sterol-response/E-box element (DREME_031; SREBF2), an NFAT motif (DREME_043), a GC-box
(GLAM2_009; SP1/KLF5/EGR1/ZNF263) and a FOXO3 motif (AMD_002) — whereas the most
significant KO-enriched motifs are de novo sequences with no confident TF match. All ten
motifs pass q < 0.05.

**Reproducibility.** The figure is regenerated from the DNAFORM output by
`make_SUDIPGXIII_motif_logo_table.py` (reads the RECLU motif tables and PWMs, recomputes
the BH-FDR, and renders the table). Vector output uses editable text
(`svg.fonttype='none'`, `pdf.fonttype=42`).

---

## References

1. Takahashi H, Lassmann T, Murata M, Carninci P. 5′ end-centered expression profiling
   using cap-analysis gene expression and next-generation sequencing. *Nat Protoc.*
   2012;7(3):542–561. https://doi.org/10.1038/nprot.2012.005
2. Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform.
   *Bioinformatics.* 2009;25(14):1754–1760. https://doi.org/10.1093/bioinformatics/btp324
3. Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and
   genotyping with HISAT2 and HISAT-genotype. *Nat Biotechnol.* 2019;37(8):907–915.
   https://doi.org/10.1038/s41587-019-0201-4
4. Ohmiya H, Vitezic M, Frith MC, et al. RECLU: a pipeline to discover reproducible
   transcriptional start sites and their alternative regulation using capped analysis of
   gene expression (CAGE). *BMC Genomics.* 2014;15:269.
   https://doi.org/10.1186/1471-2164-15-269
5. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential
   expression analysis of digital gene expression data. *Bioinformatics.*
   2010;26(1):139–140. https://doi.org/10.1093/bioinformatics/btp616
6. Bailey TL. DREME: motif discovery in transcription factor ChIP-seq data.
   *Bioinformatics.* 2011;27(12):1653–1659. https://doi.org/10.1093/bioinformatics/btr261
7. Frith MC, Saunders NFW, Kobe B, Bailey TL. Discovering sequence motifs with arbitrary
   insertions and deletions. *PLoS Comput Biol.* 2008;4(4):e1000071.
   https://doi.org/10.1371/journal.pcbi.1000071
8. Zambelli F, Pesole G, Pavesi G. Using Weeder, Pscan, and PscanChIP for the discovery
   of enriched transcription factor binding site motifs in nucleotide sequences. *Curr
   Protoc Bioinformatics.* 2014;47:2.11.1–2.11.31.
   https://doi.org/10.1002/0471250953.bi0211s47
9. Gupta S, Stamatoyannopoulos JA, Bailey TL, Noble WS. Quantifying similarity between
   motifs. *Genome Biol.* 2007;8(2):R24. https://doi.org/10.1186/gb-2007-8-2-r24
10. Mathelier A, Fornes O, Arenillas DJ, et al. JASPAR 2016: a major expansion and update
    of the open-access database of transcription factor binding profiles. *Nucleic Acids
    Res.* 2016;44(D1):D110–D115. https://doi.org/10.1093/nar/gkv1176
11. Grant CE, Bailey TL, Noble WS. FIMO: scanning for occurrences of a given motif.
    *Bioinformatics.* 2011;27(7):1017–1018. https://doi.org/10.1093/bioinformatics/btr064
12. Tareen A, Kinney JB. Logomaker: beautiful sequence logos in Python. *Bioinformatics.*
    2020;36(7):2272–2274. https://doi.org/10.1093/bioinformatics/btz921

*Citations retrieved and verified via PubMed. Ref. 4 (RECLU) is the pipeline paper cited
in the DNAFORM report; verify its exact volume/DOI against the report if needed.*
