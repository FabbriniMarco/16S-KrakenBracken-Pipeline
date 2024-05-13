#!/usr/bin/Rscript
# 10 May 2024
# Check dependancies
if (!requireNamespace("ape", quietly = TRUE))
    install.packages("ape")
if (!requireNamespace("phyloseq", quietly = TRUE))
    install.packages("phyloseq")
if (!requireNamespace("parallel", quietly = TRUE))
    install.packages("parallel")
if (!requireNamespace("doParallel", quietly = TRUE))
    install.packages("doParallel")
if (!requireNamespace("vegan", quietly = TRUE))
    install.packages("vegan")
if (!requireNamespace("otuSummary", quietly = TRUE))
    install.packages("otuSummary")
if (!requireNamespace("picante", quietly = TRUE))
    install.packages("picante")
if (!requireNamespace("openxlsx", quietly = TRUE))
    install.packages("openxlsx")
# Load libraries
library("ape")
library("phyloseq")
library("parallel")
library("doParallel")
library("otuSummary")
library("picante")
library("vegan")
library("openxlsx")
# Load data
otutab <- read.delim("final_otu_table_matrix.csv", sep=",", header=T, row.names=1)
otutab$taxa <- NULL
tree <- read.tree("/mnt/databases/kraken_db/silvaNR99/tax_slv_ssu_138.1.tre")
# Beta diversity
rooted_tree <- root(tree, 3)
unmapped <- rooted_tree$tip.label[!rooted_tree$tip.label %in% rownames(otutab)]
pruned_tree <- drop.tip(rooted_tree, unmapped)
length_tree <- compute.brlen(pruned_tree, method = "Grafen")
phydata <- phyloseq( otu_table(otutab, taxa_are_rows = TRUE), phy_tree(length_tree))
cl <- makeCluster(5)
registerDoParallel(cl)
suppressWarnings(weighted_unifrac <- UniFrac(phydata, weighted=TRUE, normalized=TRUE, parallel=TRUE))
suppressWarnings(unweighted_unifrac <- UniFrac(phydata, weighted=FALSE, normalized=TRUE, parallel=TRUE))
braycurtis <- vegdist( t(otutab), method = "bray")
jaccard <- vegdist( t(otutab), method = "jaccard")
# Alpha diversity
alphadiv = alphaDiversity(otutab, siteInCol = TRUE)
faith = pd( t(otutab), length_tree, include.root = FALSE)
colnames(faith) = c("Faith_PD", "SpeciesRichness")
alphadiv$allBio = cbind(alphadiv$allBio, faith)
#Export data
dir.create("beta_diversity")
metrics = c("braycurtis", "jaccard", "unweighted_unifrac", "weighted_unifrac")
for (m in metrics)
{
    dir.create( paste("beta_diversity", m, sep="/") )
    write.table( as.data.frame(as.matrix(get(m))), paste("beta_diversity", m, "distance-matrix.tsv", sep="/"), col.names = NA, row.names = TRUE, quote = FALSE, sep='\t' )
}
write.xlsx( alphadiv$allBio, "alpha.xlsx", colNames=TRUE, rowNames=TRUE, sheetName = "alpha_all")
wb <- loadWorkbook("alpha.xlsx")
addWorksheet(wb,"alpha_abun")
writeData(wb,"alpha_abun", alphadiv$abundBio, colNames=TRUE, rowNames=TRUE)
addWorksheet(wb,"alpha_rare") #Rare taxa are <1% relative abundance
writeData(wb,"alpha_rare", alphadiv$rareBio, colNames=TRUE, rowNames=TRUE)
