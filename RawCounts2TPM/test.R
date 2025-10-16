##############################
# NCBI项目自带文件，若没有可以根据GTF_Process自行处理得到
# NCBI的注释文件
anno <- fread("D://R//Raw2TPM//Human.GRCh38.p13.annot.tsv//Human.GRCh38.p13.annot.tsv",data.table =F)
# colnames(anno)
anno <- anno[,c("GeneID","Symbol","EnsemblGeneID","Length","GeneType")]
save(anno,file ="GRCh38.p13.annot.Rdata")
# head(anno)

############################### count2TPM
# 1.读取count
counts <- fread("D://R//Raw2TPM//GSE229904_raw_counts_GRCh38.p13_NCBI.tsv//GSE229904_raw_counts_GRCh38.p13_NCBI.tsv",data.table = F)
# counts[1:4,1:4]
rownames(counts) <- counts[,1]
counts <- counts[,-1]
# counts[1:4,1:4]


# 2.提取基因长度列
load("GRCh38.p13.annot.Rdata")
# head(anno)
rownames(anno) <- anno$GeneID
effLen <- anno[row.names(counts),"Length"]
# range(effLen)

#############################################################################
#TPM (Transcripts Per Kilobase Million) 每千个碱基的转录每百万映射读取的Transcripts
counts2TPM <-function(count=count, efflength=efflen){ 
  RPK <- count/(efflength/1000)   #每千碱基reads (reads per kilobase) 长度标准化 
  PMSC_rpk <- sum(RPK)/1e6  #RPK的每百万缩放因子 (“per million” scaling factor ) 深度标准化 
  RPK/PMSC_rpk 
}

tpm_raw <- apply(counts, 2, counts2TPM, efflength = effLen)
# tpm_raw[1:5,1:5]

# 以上都是示例代码   以下执行代码
####################################################################
path<- "D://R//yiyuzheng_prework//counts_final.rds"

counts_fin<- readRDS(path)
counts_fin<- counts_fin %>% 
  column_to_rownames("Geneid")   # 列变行名，自动删该列
counts_fin_tem <- select(counts_fin,-Length) 

tpm_raw <- apply(counts_fin_tem,2,counts2TPM ,efflength = counts_fin$Length)
# 修改列名
colnames(tpm_raw) <- sub(".*/.*-(HX[0-9]+).*", "\\1", colnames(tpm_raw))
saveRDS(tpm_raw,"tpm.rds")
