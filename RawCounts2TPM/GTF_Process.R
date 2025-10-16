#######################################################
# TGF文件处理
## 利用GenomicFeatures包导入gtf处理
txdb <- makeTxDbFromGFF("D://R//yiyuzheng_prework//Homo_sapiens.GRCh38.106.gtf//Homo_sapiens.GRCh38.106.gtf",
                        format="gtf") 
exons_gene <- exonsBy(txdb, by = "gene") ###提取基因外显子
# head(exons_gene)

##计算总外显子长度：用reduce去除掉重叠冗余的部分，,width统计长度，最后计算总长度
exons_gene_lens <- parLapply(cl,exons_gene,function(x){sum(width(reduce(x)))}) 
# exons_gene_lens[1:10]

##转换为dataframe
geneid_efflen <- data.frame(Geneid=names(exons_gene_lens),
                            Length=as.numeric(exons_gene_lens))