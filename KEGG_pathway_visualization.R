rm(list = ls())

setwd("D:/Program/RNA_SEQ/Project_2/20190219/4.Differential/1.deglist")

# org.db官方文档：http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)

rH1_H9_3vsCONTROL_up_res <- read.table("rH1_H9_3vsCONTROL/rH1_H9_3vsCONTROL_deg_up.xls",sep="\t",header=TRUE,quote="")
degs.list <- as.character(rH1_H9_3vsCONTROL_up_res$gene_name)


## id转换为EntrezID
degs.list.EntrezID <- mapIds(
  x = org.Mm.eg.db,
  keys = degs.list,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first")

# 富集
enrich.kegg.res <- enrichKEGG(gene=degs.list.EntrezID,
                              organism = "mmu",
                              keyType = "kegg")
enrich.kegg.res
barplot(enrich.kegg.res)
dotplot(enrich.kegg.res)

# ID转换
# 需要转换的Ensembl_ID
Ensembl_ID <- rH1_H9_3vsCONTROL_up_res$gene_id
# 采用bitr()函数进行转换
gene_id <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Mm.eg.db")
# 查看转换的结果
head(gene_id)

# 通路可视化
foldchanges =rH1_H9_3vsCONTROL_up_res$log2FoldChange
names(foldchanges)= gene_id$ENTREZID
head(foldchanges)

# 根据KEGG富集结果任选一条感兴趣的通路，把差异表达基因打上去
pathview(gene.data = foldchanges, pathway.id = "05164", species="mmu",kegg.native = T, same.layer = F)
pathview(gene.data = foldchanges, pathway.id = "05164", species="mmu",kegg.native = F, sign.pos="bottomleft",split.group = T)


