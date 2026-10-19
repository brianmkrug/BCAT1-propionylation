REPO_ROOT <- Sys.getenv("BCAT1_REPO", ".")  # publication: set BCAT1_REPO or run from the repository root
## Figures v2: identical to relabel_figs.R but the 4,909 Unresolved (IFN-high,
## clusters 6/18) nuclei are DROPPED ENTIRELY. 7 cell types remain, 73,151 nuclei.
## Z scores are RECOMPUTED across the retained nuclei so "Z across all nuclei" is
## honest. Palette loses the grey slot; the 7 chromatic colours already passed the
## validator --pairs all, so no revalidation needed.

suppressMessages({library(ggplot2); library(dplyr); library(ggrepel); library(patchwork)
                  library(ragg); library(hexbin); library(scales)})
set.seed(42)
SP  <- "/private/tmp/claude-501/-Users-briankrug-Desktop-BCAT1-BCAT1-scRNA-mouse/48fad6c2-80e8-4819-bcbc-9c585b6076ee/scratchpad"
RES <- file.path(REPO_ROOT, "scRNA_mouse/gsea_analysis/results")
FIG <- file.path(REPO_ROOT, "scRNA_mouse/gsea_analysis/figures")
df <- readRDS(file.path(RES, "umap_scores.rds"))

map <- c(Malignant_OPC_like="OPC", OPC="OPC", Malignant_cycling="Cycling",
         Oligodendrocyte="Oligodendrocyte", Oligodendrocyte_NFOL="Oligodendrocyte",
         Astrocyte="Astrocyte", Neuron="Neuronal",
         Microglia="Immune", Myeloid="Immune", T_cell="Immune",
         Endothelial="Vascular/stromal", VLMC_fibroblast="Vascular/stromal",
         Unresolved_IFN_high="Unresolved")
df$ct <- unname(map[df$cell_type])
df <- df[df$ct != "Unresolved", ]                       # <-- DROP unresolved entirely
cat("nuclei retained:", nrow(df), "(dropped Unresolved)\n")

LEV <- c("OPC","Cycling","Oligodendrocyte","Astrocyte","Neuronal","Immune","Vascular/stromal")
COL <- c(OPC="#009E73", Cycling="#785EF0", Oligodendrocyte="#D55E00", Astrocyte="#56B4E9",
         Neuronal="#0072B2", Immune="#CC79A7", `Vascular/stromal`="#E69F00")
LAB <- c(OPC="#008560", Cycling="#7257F0", Oligodendrocyte="#BE5400", Astrocyte="#377BA2",
         Neuronal="#0072B2", Immune="#A35D84", `Vascular/stromal`="#9A6A00")
df$ct <- factor(df$ct, levels=LEV)
SURF <- "#fcfcfb"

n <- table(df$ct)
lab_leg <- sprintf("%s (%s)", LEV, formatC(as.integer(n[LEV]), big.mark=",")); names(lab_leg) <- LEV
cat("=== cells per type (v2) ===\n"); print(n); cat("total", sum(n), "\n")

## ---------- Figure 1 v2 ----------
anchor <- df %>% group_by(ct, cell_type) %>%
  summarise(x=median(UMAP_1), y=median(UMAP_2), k=n(), .groups="drop") %>%
  group_by(ct) %>% slice_max(k, n=1) %>% ungroup() %>% mutate(label=as.character(ct))
p1 <- ggplot(df, aes(UMAP_1, UMAP_2, colour=ct)) +
  geom_point(size=0.28, alpha=0.6, shape=16) +
  scale_colour_manual(values=COL, breaks=LEV, labels=lab_leg, name=NULL,
                      guide=guide_legend(override.aes=list(size=3.2, alpha=1), ncol=1)) +
  ggrepel::geom_text_repel(data=anchor, aes(x,y,label=label), colour=LAB[as.character(anchor$ct)],
    inherit.aes=FALSE, bg.color="white", bg.r=0.18, size=4.3, fontface="bold",
    box.padding=0.8, point.padding=0.4, min.segment.length=0.25,
    segment.colour="#777777", segment.size=0.35, seed=42, max.overlaps=Inf) +
  coord_fixed() +
  labs(title="KPP/KNF mouse H3K27M glioma - cell identity",
       subtitle=paste0("73,151 nuclei | harmony-integrated UMAP | 5 samples | labelled by transcriptional cell type\n",
                       "Cells named by best-matching identity (no tumour/normal call - GFP too ambient); OPC includes formerly tumour-flagged OPC-like nuclei"),
       caption="v2: 4,909 IFN-high 'Unresolved' nuclei (clusters 6/18) dropped. Cycling = proliferative state (best lineage match OPC).") +
  theme_void(base_size=12) +
  theme(plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
        legend.position="right", legend.text=element_text(size=10, colour="#333333"),
        legend.key.height=unit(1.05,"lines"),
        plot.title=element_text(face="bold", size=15, colour="#1a1a1a"),
        plot.subtitle=element_text(size=9.5, colour="#555555", lineheight=1.25),
        plot.caption=element_text(size=8.5, colour="#777777", hjust=0), plot.margin=margin(12,12,12,12))
ggsave(file.path(FIG,"fig1_cell_identity_umap_v2.png"), p1, width=12, height=8.8, dpi=300,
       bg=SURF, device=ragg::agg_png)

## ---------- Figures 3 & 4 v2 ----------
base <- theme_void(base_size=12) +
  theme(plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
        plot.title=element_text(face="bold", size=13, colour="#1a1a1a"),
        plot.subtitle=element_text(size=9, colour="#555555"),
        legend.title=element_text(size=9), legend.text=element_text(size=8))

make_z <- function(col, nice, outfile, ngenes, pct0) {
  df$Z <- as.numeric(scale(df[[col]]))                  # recomputed on retained nuclei
  p_z <- ggplot(df, aes(UMAP_1, UMAP_2, z=Z)) +
    stat_summary_hex(bins=90, fun=mean) +
    scale_fill_gradient2(low="#2166AC", mid="#F0F0F0", high="#B2182B", midpoint=0,
                         limits=c(-1,1), oob=scales::squish, name="mean Z") +
    coord_fixed() +
    labs(title=paste0(nice, " - Z score of per-nucleus UCell"),
         subtitle=paste0("diverging ramp, neutral grey midpoint | hex-binned mean, clipped to +/-1 SD | ",
                         pct0, "% score 0 (", ngenes, " genes) | Unresolved dropped")) +
    base + theme(legend.position="right")
  ord <- df %>% group_by(ct) %>% summarise(m=mean(Z), .groups="drop") %>% arrange(desc(m))
  d2 <- df; d2$cto <- factor(d2$ct, levels=ord$ct)
  st <- d2 %>% group_by(cto) %>% summarise(n=n(), mean_Z=mean(Z), median_Z=median(Z), .groups="drop")
  p_box <- ggplot(d2, aes(cto, Z, fill=cto)) +
    geom_hline(yintercept=0, colour="#cccccc", linewidth=0.4) +
    geom_boxplot(outlier.shape=NA, linewidth=0.3, colour="#444444", width=0.66) +
    geom_point(data=st, aes(cto, mean_Z), inherit.aes=FALSE, shape=23, size=2.3,
               fill="white", colour="#111111", stroke=0.6) +
    scale_fill_manual(values=COL, guide="none") +
    scale_x_discrete(labels=function(x) paste0(x, "\n(n=", formatC(st$n[match(x,st$cto)], big.mark=","), ")")) +
    coord_cartesian(ylim=c(-2, 2.6)) +
    labs(title=paste0(nice, " score by cell type"),
         subtitle="per-nucleus UCell, Z-scored across 73,151 retained nuclei | box = median/IQR | white diamond = mean | ordered by mean",
         x=NULL, y="Z score") +
    theme_minimal(base_size=12) +
    theme(plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          panel.grid.major.y=element_line(colour="#e8e8e8", linewidth=0.3),
          axis.text.x=element_text(size=8.5, colour="#333333", lineheight=1.1),
          plot.title=element_text(face="bold", size=13), plot.subtitle=element_text(size=8.5, colour="#555555"))
  fig <- p_z / p_box + plot_layout(heights=c(1.45,1)) +
    plot_annotation(theme=theme(plot.background=element_rect(fill=SURF,colour=NA)))
  ggsave(file.path(FIG,outfile), fig, width=12, height=11, dpi=300, bg=SURF, device=ragg::agg_png)
  cat("\n===", nice, "(v2) Z by cell type ===\n")
  print(as.data.frame(st %>% mutate(across(c(mean_Z,median_Z), ~round(.x,3)))), row.names=FALSE)
}
make_z("Cholesterol_synthesis", "Cholesterol synthesis", "fig3_cholesterol_zscore_by_celltype_v2.png", 21, 35.7)
make_z("Fatty_acid_metabolism", "Fatty acid metabolism", "fig4_fattyacid_zscore_by_celltype_v2.png", 152, 0.2)
cat("\nDONE\n")
