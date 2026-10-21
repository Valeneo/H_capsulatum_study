# DAPC (Discriminant Analysis of Principal Components) on the LD-pruned SNP set.
# Run: Rscript scripts/06b_dapc.R
suppressPackageStartupMessages({library(vcfR); library(adegenet); library(ggplot2)})

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
out  <- file.path(root, "results/06_popgen")

# The LD-pruned SNPs (same set as PCA/ADMIXTURE), written by 06_popgen.sh, genotypes still haploid.
vcf  <- read.vcfR(file.path(out, "pruned.haploid.vcf.gz"), verbose = FALSE)

# haploid organism -> ploidy = 1
gl <- vcfR2genlight(vcf); ploidy(gl) <- 1

# 1. find.clusters: k-means on the principal components for k = 1..10, pick k by BIC (lowest / elbow).
#    This is UNSUPERVISED - no need to say in advance which isolates belong together.
set.seed(1)
fc <- find.clusters(gl, max.n.clust = 10, n.pca = 100, choose.n.clust = FALSE, criterion = "diffNgroup")
grp <- fc$grp
cat("clusters chosen:", length(unique(grp)), "\n"); print(table(grp))

# 2. DAPC: finds the axes that BEST SEPARATE the groups (PCA only finds axes of most variance).
#    Too many PCs overfits; xval.dapc picks the number that predicts membership best on held-out isolates.
xv <- xvalDapc(tab(gl, NA.method = "mean"), grp, n.pca.max = 100, training.set = 0.9, result = "overall", n.rep = 20, parallel = "no")
dp <- dapc(gl, grp, n.pca = as.integer(xv[["Number of PCs Achieving Lowest MSE"]]), n.da = length(unique(grp)) - 1)

write.table(data.frame(run = indNames(gl), dapc_group = grp, dp$ind.coord), file.path(out, "dapc_assignments.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
png(file.path(out, "dapc_scatter.png"), 1600, 1200, res = 200); scatter(dp, scree.da = TRUE, legend = TRUE); dev.off()
png(file.path(out, "dapc_membership.png"), 2000, 800, res = 200); compoplot(dp, legend = TRUE); dev.off()
