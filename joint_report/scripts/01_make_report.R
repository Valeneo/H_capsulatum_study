# Builds tables and composite figures for the clade-wise SNP results.
# Inputs (read-only): <clade>_snpeff_summary.csv, <clade>_snpeff_summary.genes.txt, <clade>_joint.ann.vcf.gz
#                     report/tables/raw/<clade>_bcftools_stats.txt  (made by 00_bcftools_stats.sh)
# Outputs: report/tables/*.csv, report/figures/*.png|pdf
# Run: Rscript report/scripts/01_make_report.R
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(stringr); library(forcats); library(patchwork); library(scales); library(ggrepel)
})

args <- commandArgs(FALSE)
sf   <- sub("--file=", "", args[grep("--file=", args)])
ROOT <- normalizePath(file.path(dirname(sf), "..", ".."))
FIG  <- file.path(ROOT, "report", "figures"); TAB <- file.path(ROOT, "report", "tables")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE); dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
BCF  <- Sys.which("bcftools"); if (BCF == "") stop("bcftools not found on PATH")

clades <- c("NAm1", "African", "LatinAm", "G217B")
clade_label <- c(NAm1 = "NAm1 (ref WU24)", African = "African", LatinAm = "Latin American", G217B = "G217B")
pal <- c(NAm1 = "#0072B2", African = "#E69F00", LatinAm = "#009E73", G217B = "#CC79A7")   # Okabe-Ito
impact_cols <- c(HIGH = "#D55E00", MODERATE = "#E69F00", LOW = "#56B4E9", MODIFIER = "#999999")

theme_set(theme_classic(base_size = 11) +
  theme(plot.tag = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = "grey35", size = 9), legend.key.size = unit(0.4, "cm"),
        axis.line = element_line(colour = "grey30"), strip.background = element_blank(),
        strip.text = element_text(face = "bold")))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, device = dev)
}

## ---------- parsers -------------------------------------------------------------------------
section <- function(lines, title) {
  i <- which(trimws(lines) == paste0("# ", title)); if (!length(i)) return(NULL)
  i <- i[1] + 1; while (i <= length(lines) && trimws(lines[i]) == "") i <- i + 1
  j <- i; while (j <= length(lines) && trimws(lines[j]) != "" && !startsWith(lines[j], "#")) j <- j + 1
  if (j <= i) NULL else lines[i:(j - 1)]
}
tab_section <- function(lines, title) {
  s <- section(lines, title); if (is.null(s) || length(s) < 2) return(tibble())
  parts <- lapply(strsplit(s, ",", fixed = TRUE), trimws); hdr <- parts[[1]]
  m <- do.call(rbind, lapply(parts[-1], function(p) p[seq_along(hdr)])); df <- as_tibble(m, .name_repair = "minimal"); names(df) <- hdr; df
}
kv <- function(lines, key) {
  l <- lines[startsWith(trimws(lines), key)]; if (!length(l)) return(NA_real_)
  as.numeric(gsub("[^0-9.]", "", strsplit(sub("^[^,]*,", "", l[1]), ",")[[1]][1]))
}
vcf_header <- function(clade) system2(BCF, c("view", "-h", file.path(ROOT, paste0(clade, "_joint.ann.vcf.gz"))), stdout = TRUE)

## ---------- per-clade data ---------------------------------------------------------------------
info <- list(); impact <- list(); fclass <- list(); region <- list(); effects <- list(); genes <- list(); psc <- list(); afs <- list(); tstv <- list()
for (c in clades) {
  L <- readLines(file.path(ROOT, paste0(c, "_snpeff_summary.csv")), warn = FALSE)
  H <- vcf_header(c); ctg <- H[startsWith(H, "##contig")]
  ref_len <- sum(as.numeric(sub(".*length=([0-9]+).*", "\\1", ctg)))
  ref_asm <- sub(".*assembly=([^,>]+).*", "\\1", ctg[1])
  db <- trimws(sub("^[^,]*,", "", L[startsWith(L, "Genome ,")][1]))
  stats_f <- file.path(TAB, "raw", paste0(c, "_bcftools_stats.txt")); S <- readLines(stats_f, warn = FALSE)

  p <- read_tsv(I(S[startsWith(S, "PSC\t")]), col_names = c("tag","id","sample","nRefHom","nNonRefHom","nHets","nTs","nTv","nIndels","avgDP","nSingletons","nHapRef","nHapAlt","nMissing"), show_col_types = FALSE)
  psc[[c]] <- p |> transmute(clade = c, sample, mean_depth = avgDP, hap_ref = nHapRef, hap_alt = nHapAlt, missing = nMissing,
                             pct_nonref = 100 * nHapAlt / (nHapRef + nHapAlt), pct_missing = 100 * nMissing / (nHapRef + nHapAlt + nMissing))
  a <- read_tsv(I(S[startsWith(S, "AF\t")]), col_names = c("tag","id","af","n_snps","ts","tv","n_indels","rc","ri","na"), show_col_types = FALSE)
  afs[[c]] <- a |> transmute(clade = c, af, n_snps)
  ts <- read_tsv(I(S[startsWith(S, "TSTV\t")]), col_names = c("tag","id","ts","tv","tstv","ts1","tv1","tstv1"), show_col_types = FALSE)
  tstv[[c]] <- tibble(clade = c, ts_tv = ts$tstv[1])

  variants <- kv(L, "Number_of_variants_processed"); n_eff <- kv(L, "Number_of_effects")
  err_all <- identical(n_eff, variants) && is.na(kv(L, "Missense_Silent_ratio")) | nrow(tab_section(L, "Effects by impact")) == 0
  info[[c]] <- tibble(clade = c, snpeff_db = db, reference_assembly = ref_asm, reference_contigs = length(ctg), reference_length_bp = ref_len,
                      isolates = nrow(p), snps = variants, snps_per_kb = variants / (ref_len / 1000), snps_per_isolate_median = median(p$nHapAlt),
                      effects = n_eff, warnings = kv(L, "Warnings"), warnings_pct_of_snps = 100 * kv(L, "Warnings") / variants,
                      missense_silent = kv(L, "Missense_Silent_ratio"), ts_tv = ts$tstv[1], annotation_ok = !err_all)
  strip_pct <- function(df) df |> mutate(Count = as.numeric(Count), Percent = as.numeric(gsub("%", "", Percent)))
  t1 <- tab_section(L, "Effects by impact");            if (nrow(t1)) impact[[c]]  <- strip_pct(t1) |> mutate(clade = c)
  t2 <- tab_section(L, "Effects by functional class");  if (nrow(t2)) fclass[[c]]  <- strip_pct(t2) |> mutate(clade = c)
  t3 <- tab_section(L, "Count by genomic region");      if (nrow(t3)) region[[c]]  <- strip_pct(t3) |> mutate(clade = c)
  t4 <- tab_section(L, "Count by effects");             if (nrow(t4)) effects[[c]] <- strip_pct(t4) |> mutate(clade = c)

  gf <- file.path(ROOT, paste0(c, "_snpeff_summary.genes.txt")); G <- readLines(gf, warn = FALSE)
  hi <- grep("^#GeneName", G); body <- G[!startsWith(G, "#")]
  if (length(hi) && length(body)) genes[[c]] <- read_tsv(I(c(sub("^#", "", G[hi]), body)), show_col_types = FALSE, progress = FALSE) |> mutate(clade = c)
}
info <- bind_rows(info) |> mutate(clade = factor(clade, clades)); impact <- bind_rows(impact); fclass <- bind_rows(fclass)
region <- bind_rows(region); effects <- bind_rows(effects); genes <- bind_rows(genes); psc <- bind_rows(psc); afs <- bind_rows(afs); tstv <- bind_rows(tstv)

write_csv(info, file.path(TAB, "clade_summary.csv")); write_csv(impact, file.path(TAB, "effects_by_impact.csv"))
write_csv(fclass, file.path(TAB, "effects_by_functional_class.csv")); write_csv(region, file.path(TAB, "effects_by_region.csv"))
write_csv(effects, file.path(TAB, "effects_by_type.csv")); write_csv(psc, file.path(TAB, "per_isolate_summary.csv"))
write_csv(tstv, file.path(TAB, "ts_tv_by_clade.csv"))
groups <- psc |> mutate(divergence_group = if_else(pct_nonref < 5, "near-identical to reference (<5% non-ref)", "divergent (>=5% non-ref)")) |> count(clade, divergence_group)
write_csv(groups, file.path(TAB, "divergence_groups_by_clade.csv"))

excluded <- tribble(~Run_Accession, ~country, ~mapping_rate_reported, ~reference_in_that_run, ~note,
  "SRR20306115", "China",        "55.7%", "WU24", "absent from all four clade sets; reason not recorded in the results folder",
  "SRR29905426", "South Africa", "65.6%", "H88",  "absent from all four clade sets; reason not recorded in the results folder",
  "SRR29905402", "Egypt",        "67.2%", "H88",  "absent from all four clade sets; reason not recorded in the results folder")
write_csv(excluded, file.path(TAB, "isolates_not_in_any_clade_set.csv"))

top_high <- genes |> arrange(clade, desc(variants_impact_HIGH)) |> group_by(clade) |> slice_head(n = 15) |>
  transmute(clade, GeneName, GeneId, high = variants_impact_HIGH, moderate = variants_impact_MODERATE, low = variants_impact_LOW,
            stop_gained = variants_effect_stop_gained, start_lost = variants_effect_start_lost, splice_donor = variants_effect_splice_donor_variant, splice_acceptor = variants_effect_splice_acceptor_variant)
write_csv(top_high, file.path(TAB, "top15_high_impact_genes_per_clade.csv"))
gene_totals <- genes |> group_by(clade) |> summarise(genes_in_table = n(), genes_with_HIGH = sum(variants_impact_HIGH > 0), genes_with_MODERATE = sum(variants_impact_MODERATE > 0), .groups = "drop")
write_csv(gene_totals, file.path(TAB, "gene_level_totals.csv"))

ann_clades <- info$clade[info$annotation_ok] |> as.character()
cl <- function(x) factor(clade_label[as.character(x)], levels = clade_label)
fill_cl <- scale_fill_manual(values = setNames(pal, clade_label), guide = "none")
col_cl  <- scale_colour_manual(values = setNames(pal, clade_label), name = NULL)

## ---------- Figure 1: overview ---------------------------------------------------------------
pA <- ggplot(info, aes(cl(clade), isolates, fill = cl(clade))) + geom_col(width = 0.65) +
  geom_text(aes(label = isolates), vjust = -0.4, size = 3.4) + fill_cl +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Isolates", title = "Isolates per clade set", subtitle = paste0(sum(info$isolates), " of 183 isolates; 3 are in no set"))
ref_kb <- setNames(info$reference_length_bp / 1000, as.character(info$clade))
pB <- psc |> mutate(per_kb = hap_alt / ref_kb[clade], clade_l = cl(clade)) |>
  ggplot(aes(clade_l, per_kb, colour = clade_l)) + geom_jitter(width = 0.18, size = 1.8, alpha = 0.8) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, middle.linewidth = 1.5, colour = "grey20", fill = NA) + col_cl +
  labs(x = NULL, y = "Non-reference SNPs per kb of reference", title = "Divergence from the clade reference",
       subtitle = "one dot per isolate, bar = median; sets with two clusters mix lineages") + theme(legend.position = "none")
imp <- impact |> filter(Type %in% names(impact_cols)) |> mutate(Type = factor(Type, names(impact_cols)), clade_l = cl(clade))
pC <- ggplot(imp, aes(Percent, Type, colour = clade_l)) + geom_point(size = 3, position = position_dodge(width = 0.55)) +
  scale_x_log10(labels = label_number(drop0trailing = TRUE), breaks = c(0.01, 0.1, 1, 10, 100)) + col_cl +
  labs(x = "% of predicted effects (log scale)", y = NULL, title = "Predicted impact of SNP effects",
       subtitle = "Latin American not annotated (see Figure 5)") + theme(legend.position = "bottom")
pD <- fclass |> filter(Type %in% c("MISSENSE", "SILENT", "NONSENSE")) |> mutate(clade_l = cl(clade), Type = str_to_title(Type)) |>
  filter(Type != "Nonsense") |> select(clade_l, Type, Count) |>
  ggplot(aes(clade_l, Count / 1000, fill = Type)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(Missense = "#D55E00", Silent = "#56B4E9"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_comma()) +
  labs(x = NULL, y = "Coding SNP effects (thousands)", title = "Missense versus silent",
       subtitle = paste0("missense:silent = ", paste0(clade_label[as.character(info$clade[info$annotation_ok])], " ", number(info$missense_silent[info$annotation_ok], 0.01), collapse = "; "))) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 8, colour = "grey35"))
fig1 <- (pA | pB) / (pC | pD) + plot_annotation(tag_levels = "A", title = "Clade-wise SNP results against clade-specific references")
save_fig(fig1, "Figure1_overview", 11, 9)

## ---------- Figure 2: genomic context ----------------------------------------------------------
reg <- region |> mutate(Type = case_when(str_detect(Type, "SPLICE") ~ "SPLICE", str_detect(Type, "UTR") ~ "UTR", TRUE ~ Type)) |> group_by(clade, Type) |> summarise(Count = sum(Count), .groups = "drop") |>
  group_by(clade) |> mutate(pct = 100 * Count / sum(Count)) |> ungroup() |>
  mutate(Type = factor(str_to_title(Type), levels = c("Upstream", "Downstream", "Intergenic", "UTR", "Exon", "Intron", "Splice")), clade_l = cl(clade))
pA2 <- ggplot(reg, aes(pct, fct_rev(clade_l), fill = Type)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c(Upstream = "#BBBBBB", Downstream = "#888888", Intergenic = "#DDDDDD", UTR = "#009E73", Exon = "#0072B2", Intron = "#56B4E9", Splice = "#D55E00"), name = NULL) +
  scale_x_continuous(expand = c(0, 0)) + labs(x = "% of predicted effects", y = NULL, title = "Where SNPs fall in the genome", subtitle = "NAm1 annotation has no UTR class, so its intergenic share is not comparable")
ef_keep <- c("missense_variant", "synonymous_variant", "stop_gained", "start_lost", "stop_lost", "splice_donor_variant", "splice_acceptor_variant")
ef <- effects |> filter(Type %in% ef_keep) |> left_join(select(info, clade, snps), by = "clade") |>
  mutate(per10k = 1e4 * Count / snps, Type = fct_reorder(str_replace_all(Type, "_", " "), per10k, .desc = TRUE), clade_l = cl(clade))
pB2 <- ggplot(ef, aes(per10k, fct_rev(Type), colour = clade_l)) + geom_point(size = 3, position = position_dodge(width = 0.55)) +
  scale_x_log10(labels = label_comma()) + col_cl + labs(x = "Events per 10,000 SNP sites (log scale)", y = NULL, title = "Coding-region effects", subtitle = "Latin American not annotated") +
  theme(legend.position = "bottom")
fig2 <- pA2 / pB2 + plot_layout(heights = c(1, 1.5)) + plot_annotation(tag_levels = "A")
save_fig(fig2, "Figure2_genomic_context", 8, 8)

## ---------- Figure 3: genes with the most HIGH-impact variants ---------------------------------
gene_plot <- function(c) {
  d <- genes |> filter(clade == c) |> arrange(desc(variants_impact_HIGH)) |> slice_head(n = 12) |>
    mutate(label = if_else(is.na(GeneName) | GeneName == "" | GeneName == GeneId, str_remove(GeneId, "^gene-"), GeneName))
  d$label <- make.unique(d$label)
  ggplot(d, aes(variants_impact_HIGH, fct_reorder(label, variants_impact_HIGH))) + geom_col(fill = pal[[c]], width = 0.7) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08)), breaks = function(x) unique(floor(pretty(x)))) +
    labs(x = "HIGH-impact variant sites", y = NULL, title = clade_label[[c]], subtitle = paste0(sum(genes$variants_impact_HIGH[genes$clade == c] > 0), " genes with at least one"))
}
fig3 <- wrap_plots(lapply(intersect(clades, unique(genes$clade)), gene_plot), nrow = 1) +
  plot_annotation(tag_levels = "A", title = "Genes carrying the most HIGH-impact variants",
                  subtitle = "SnpEff flagged many gene models as incomplete or with internal stop codons (Figure 5), so treat these as candidates only")
save_fig(fig3, "Figure3_top_high_impact_genes", 13, 5)

## ---------- Figure 4: per-isolate variation --------------------------------------------------
psc <- psc |> mutate(clade_l = cl(clade))
out_lab <- function(d, v) d |> group_by(clade_l) |> mutate(is_out = {s <- boxplot.stats(.data[[v]])$out; .data[[v]] %in% s}) |> ungroup() |> filter(is_out)
pts <- function(y, ylab, ttl, logy = FALSE, lab = NULL) {
  p <- ggplot(psc, aes(clade_l, .data[[y]], colour = clade_l)) +
    geom_jitter(width = 0.18, size = 1.8, alpha = 0.8) +
    stat_summary(fun = median, geom = "crossbar", width = 0.55, middle.linewidth = 1.5, colour = "grey20", fill = NA) +
    col_cl + labs(x = NULL, y = ylab, title = ttl) + theme(legend.position = "none")
  if (logy) p <- p + scale_y_log10()
  if (!is.null(lab)) p <- p + geom_text_repel(data = filter(psc, .data[[y]] > lab), aes(label = sample), size = 2.3, max.overlaps = 30, segment.size = 0.2, colour = "grey20")
  p
}
pA4 <- pts("pct_nonref", "Non-reference calls (% of called sites)", "Divergence from the clade reference")
pB4 <- pts("pct_missing", "Missing genotypes (% of sites)", "Missing data (labelled above 15%)", lab = 15)
pC4 <- pts("mean_depth", "Mean depth at SNP sites (log scale)", "Sequencing depth", logy = TRUE)
fig4 <- (pA4 | pB4 | pC4) + plot_annotation(tag_levels = "A", title = "Variation among isolates within each clade set", subtitle = "dots are isolates; the bar is the median. Two-cluster patterns mean a set contains more than one lineage")
save_fig(fig4, "Figure4_per_isolate", 14, 5)

## ---------- Figure 5: reliability of the references and annotation ---------------------------
pA5 <- ggplot(info, aes(cl(clade), reference_contigs, fill = cl(clade))) + geom_col(width = 0.65) + fill_cl +
  geom_text(aes(label = comma(reference_contigs)), vjust = -0.4, size = 3.2) +
  scale_y_log10(expand = expansion(mult = c(0, 0.25)), labels = label_comma()) +
  labs(x = NULL, y = "Contigs in the reference (log scale)", title = "Reference fragmentation", subtitle = "NAm1 uses the complete WU24 genome")
w5 <- info |> mutate(status = if_else(annotation_ok, "Warnings", "Annotation failed"), val = if_else(annotation_ok, warnings_pct_of_snps, 100), clade_l = cl(clade))
pB5 <- ggplot(w5, aes(clade_l, val, fill = status)) + geom_col(width = 0.65) +
  geom_text(aes(label = paste0(number(val, 0.1), "%")), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(Warnings = "#E69F00", `Annotation failed` = "#D55E00"), name = NULL) +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
  labs(x = NULL, y = "% of SNPs", title = "Annotation reliability", subtitle = "gene-model warnings; Latin American: no contig matched the database") + theme(legend.position = "bottom")
pC5 <- ggplot(tstv |> mutate(clade_l = cl(clade)), aes(clade_l, ts_tv, fill = clade_l)) + geom_col(width = 0.65) + fill_cl +
  geom_text(aes(label = number(ts_tv, 0.01)), vjust = -0.4, size = 3.2) + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Transition / transversion ratio", title = "Ts/Tv of SNP sites")
sfs <- afs |> group_by(clade) |> mutate(frac = n_snps / sum(n_snps), bin = cut(af, breaks = c(-Inf, 0.1, 0.2, 0.4, 0.6, 0.8, 1), labels = c("<=10%", "10-20%", "20-40%", "40-60%", "60-80%", ">80%"))) |> ungroup() |>
  group_by(clade, bin) |> summarise(frac = sum(frac), .groups = "drop") |> filter(frac > 0) |> mutate(clade_l = cl(clade))
pD5 <- ggplot(sfs, aes(bin, 100 * frac, colour = clade_l, group = clade_l)) + geom_line(linewidth = 0.7) + geom_point(size = 2) + col_cl +
  scale_y_log10(labels = label_number(drop0trailing = TRUE, suffix = "%"), breaks = c(0.01, 0.1, 1, 10, 100)) +
  labs(x = "Isolates carrying the non-reference allele (from INFO/AF)", y = "% of SNP sites (log scale)", title = "Allele frequency spectrum",
       subtitle = "most sites are rare; empty bins are omitted") + theme(legend.position = "bottom")
fig5 <- (pA5 | pB5) / (pC5 | pD5) + plot_annotation(tag_levels = "A", title = "Quality of references and annotation")
save_fig(fig5, "Figure5_reliability", 11, 9)

## ---------- console summary ----------------------------------------------------------------------
print(as.data.frame(info |> select(clade, isolates, snps, snps_per_kb, reference_contigs, warnings_pct_of_snps, missense_silent, ts_tv, annotation_ok)), digits = 4)
cat("\nfigures:", paste(list.files(FIG, "\\.png$"), collapse = ", "), "\n")
