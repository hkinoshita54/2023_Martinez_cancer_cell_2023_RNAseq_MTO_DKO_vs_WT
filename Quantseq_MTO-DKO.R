## Load packages ----
library(tximport)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)   # for mouse cells
library(tidyverse)
library(fgsea)

## Make sample_table ----
sample_table = read_csv("SraRunTable.txt") %>%  ## read_csv returns tibble
  select(`Sample Name`, Genotype) %>%
  unique()
# View(sample_table)

## make objects for tximport arguments
sample_files = paste0(pull(sample_table, `Sample Name`), "/quant.sf")
names(sample_files) = pull(sample_table, `Sample Name`)

## gene_map.csv is created in command line from the gencode file
gene_map_mm = read_csv("gene_map_mm.csv", col_names = c("ensmustid", "ensmusgid"))
# View(gene_map_mm)

## make count data WITHOUT NORMALIZATION
count_data = tximport(files = sample_files,
                      type = "salmon",
                      tx2gene = gene_map_mm,
                      ignoreTxVersion = TRUE)

## DESeqDataSetFromTximort ----
## modify sample_table for colData argument of DESeqDataSetFromTximort
sample_table = as.data.frame(sample_table)  ## DESeq2 doesn't recognize tibble
colnames(sample_table)[1] = "sample_name"  ## to remove space from the column name

## set proper conditions for the experiment
conditions = c(rep("DKO", times = 3), rep("WT", times = 3))
conditions = factor(conditions, levels = c("DKO", "WT"))
sample_table$conditions = conditions

deseq_dataset = DESeqDataSetFromTximport(txi = count_data,
                                         colData = sample_table,
                                         design = ~conditions)

# take a look at count data. DESeq2 rounds count to integers
# counts(deseq_dataset)[1:6, 1:3] ; count_data$counts[1:6, 1:3]

## DESeq() ---
## includes estimateSizeFactors(), estimateDispersions() and nbinomWaldTest()
deseq_dataset = DESeq(deseq_dataset)

## PCA ----
vst = varianceStabilizingTransformation(deseq_dataset)
plotPCA(vst, intgroup = "conditions")

## Retrieve and filter results ----
results = results(deseq_dataset, contrast = c("conditions", "DKO", "WT"))
results_df = data.frame(results)
results_df = rownames_to_column(results_df, var = "ensmusg")
results_filter1 = filter(results_df, complete.cases(results))

results_filter2 = filter(results_filter1, padj < 0.05)
results_filter3 = filter(results_filter2, abs(log2FoldChange) > 1)

## Prepare top 500 up- and down-regulated genes for lisa cistrome ----
DKO_up_500 <- results_filter3 %>% 
  filter(log2FoldChange > 0) %>% 
  arrange(padj) %>% 
  slice_head(n = 500)
write.csv(DKO_up_500, file = "DKO_up_500.csv")

DKO_dn_500 <- results_filter3 %>% 
  filter(log2FoldChange < 0) %>% 
  arrange(padj) %>% 
  slice_head(n = 500)
write.csv(DKO_dn_500, file = "DKO_dn_500.csv")

## MA plot ----
plotMA(results)



# annotation ----
## use proper version of ensembl for the gencode transcriptome
ensembl110 = useEnsembl(biomart = "ensembl", version = 110)
# View(listDatasets(ensembl110))
ensembl110 = useDataset("mmusculus_gene_ensembl", mart = ensembl110)
# View(listAttributes(ensembl110))
# View(listFilters(ensembl110))
annotation_mm = getBM(attributes = c("ensembl_gene_id",
                                     "chromosome_name",
                                     "start_position",
                                     "end_position",
                                     "strand",
                                     "gene_biotype",
                                     "external_gene_name",
                                     "description"),
                      filters = c("ensembl_gene_id"),
                      values = results_filter1$ensmusg,
                      mart = ensembl110)

## connect annotation to gene matrix
results_anno_df = left_join(results_filter1, annotation_mm, by = c("ensmusg" = "ensembl_gene_id"))
View(results_anno_df)
write.csv(results_anno_df, file = "results_anno_df.csv")

## volcano plot ----
## set threshold of logFC and padj
results_anno_df$test = results_anno_df$padj < 0.05 & abs(results_anno_df$log2FoldChange) > 1

g = ggplot(results_anno_df, aes(x = log2FoldChange, y = -log10(padj), name = external_gene_name)) +
  geom_point(aes(color = test), size = 1, alpha = 0.3) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = 1, color = "green", linetype = 3) +
  geom_vline(xintercept = -1, color = "green", linetype = 3) +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = 3) +
  # xlim(-3, 3) +
  ylim(0, 100) +
  theme_bw() +
  theme(legend.position = "none")
g
# interactive viewer
# ggplotly(g)

## export gene matrix for GSEA outside of R ----
## ENSMUSG was not replaced with gene name
## low count genes must be filtered out
gene_matrix = counts(deseq_dataset, normalized = TRUE) %>%
  data.frame %>%
  rownames_to_column(var = "ensmusg") %>%
  filter(complete.cases(results))
View(gene_matrix)
write.csv(gene_matrix, file = "gene_matrix.csv")

## annotated (Ensembl id replaced with external gene name)
annotation_mm_gene_name <- annotation_mm %>% select(ensembl_gene_id, external_gene_name)
gene_matrix_annotated <- left_join(gene_matrix, annotation_mm_gene_name, by = c("ensmusg" = "ensembl_gene_id")) %>% 
  select(-ensmusg) %>%
  select(external_gene_name, everything()) %>% 
  filter(!(external_gene_name == "")) %>% 
  distinct(external_gene_name, .keep_all = TRUE)
write.csv(gene_matrix_annotated, file = "gene_matrix_annotated.csv")

## Heatmap ----
results_anno2 = filter(results_anno_df, padj < 0.05)
results_anno3 = filter(results_anno2, abs(log2FoldChange) > 1)
degs = results_anno3$ensmusg
vst_mat = assay(vst)
results_hm = vst_mat[degs,]
rownames(results_hm) = results_anno3$external_gene_name
heatmap(results_hm)
pheatmap(results_hm, fontsize_row = 4, scale = "row")

## GO enrichment ----
entrez = getBM(attributes = c("entrezgene_id"),
               filters = c("ensembl_gene_id"),
               values = results_anno3$ensmusg,
               mart = ensembl110)
entrez = entrez$entrezgene_id
entrez = as.character(entrez)

universe = getBM(attributes = c("entrezgene_id"),
                 filters = c("ensembl_gene_id"),
                 values = results_anno_df$ensmusg,
                 mart = ensembl110)
universe = universe$entrezgene_id
universe = as.character(universe)

ego = enrichGO(gene = entrez,
               OrgDb = org.Mm.eg.db,
               ont = "BP",
               universe = universe,
               readable = TRUE)

barplot(ego, showCategory = 20)
dotplot(ego, showCategory = 20)

fold_changes = results_anno3$log2FoldChange
names(fold_changes) = results_anno3$external_gene_name
cnetplot(ego,
         showCategory = 10,
         foldChange = fold_changes)
goplot(ego)


## fgsea ----

# create a table to map mouse gene IDs to human symbol
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() 
bm[bm == ""] <- NA
bm <- bm %>% na.omit()

# join to the result table
res <- inner_join(results_anno_df, bm, by = c("ensmusg"="ensembl_gene_id"))

# select only test statistic and gene symbol
# remove NA
# collapse genes (in case you have multiple test statistics for the same symbol)
res2 <- res %>% 
  dplyr::select(hsapiens_homolog_associated_gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  summarize(stat=mean(stat))

# make ranked gene list by using tibble::deframe
ranks <- deframe(res2)

# first you need to download .gmt file from MSigDB
# Load the pathways into a named list
pathways.hallmark <- gmtPathways("./h.all.v2023.1.Hs.symbols.gmt")

# run fgsea with 1000 permutations
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

# tidy the results
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# plot NES. color shows if padj < 0.25
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

# enrichment plot
plotEnrichment(pathways.hallmark[["HALLMARK_TGF_BETA_SIGNALING"]],
               ranks) + labs(title="HALLMARK_TGF_BETA_SIGNALING")

# table plot
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.hallmark[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)