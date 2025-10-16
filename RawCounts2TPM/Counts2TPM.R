#############################################
# step1 运行GTF_Prcess，得到geneid_efflen文件
geneid_efflen[1:10]

##########################################
# step2 读取counts文件并进行相应处理 根据counts数据类型使用不同的读取方法
path<- "D://R//yiyuzheng_prework//counts_final.rds"

counts_mer<- readRDS(path)
counts_mer<- counts_fin %>% 
  column_to_rownames("Geneid")   # 列变行名，自动删该列

#########################################
# step3 将counts文件和geneid_efflen文件进行ID匹配

# 根据tgf处理得到的基因长度进行匹配
counts_fin <- merge(counts_mer, geneid_efflen, by = "Geneid", all.x = TRUE)

#########################################
# step4 运行counts2TPM函数

#TPM (Transcripts Per Kilobase Million) 每千个碱基的转录每百万映射读取的Transcripts
# 函数本体
counts2TPM <-function(count=count, efflength=efflen){ 
  RPK <- count/(efflength/1000)   #每千碱基reads (reads per kilobase) 长度标准化 
  PMSC_rpk <- sum(RPK)/1e6  #RPK的每百万缩放因子 (“per million” scaling factor ) 深度标准化 
  RPK/PMSC_rpk 
}
# 函数运行
counts_fin_tem <- select(counts_fin,-Length) #选定特定行 减小运行内存
tpm_raw <- apply(counts_fin_tem,2,counts2TPM ,efflength = counts_fin$Length)

#############################################
# step5 微调格式并保存
# 修改列名
colnames(tpm_raw) <- sub(".*/.*-(HX[0-9]+).*", "\\1", colnames(tpm_raw))
saveRDS(tpm_raw,"tpm.rds")