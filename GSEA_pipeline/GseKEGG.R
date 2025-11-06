
# 构建genelist
# 假设 data_file 列名：Gene（SYMBOL），logFC
map_df <- bitr(data_file$Gene,
               fromType = "SYMBOL",
               toType   = "ENTREZID",
               OrgDb    = org.Hs.eg.db)

# 把原始 logFC 合并进来，并按 EntrezID 排好序
geneList_ENTREZ <- data_file %>%
  inner_join(map_df, by = c("Gene" = "SYMBOL")) %>%   # 丢没 map 上的
  select(ENTREZID, logFC) %>%                         # 只要这两列
  deframe() %>%                                       # → named numeric
  sort(decreasing = TRUE)                             # GSEA 必须降序


#########################################################################
##KEGG 富集分析（gseKEGG）------------------------------------------
kegg_gsea <- gseKEGG(
  geneList      = geneList_ENTREZ,      # ① 已排好序的 named numeric 向量；名字必须是 ENTREZID！
  organism      = "hsa",         # ② 人类通路前缀；小鼠用 "mmu"，大鼠 "rno"
  minGSSize     = 10,            # ③ 通路至少 10 个基因，太小没统计意义
  maxGSSize     = 500,           # ④ 超过 500 个基因的通路不要（太泛，稀释信号）
  pvalueCutoff  = 0.05,          # ⑤ 原始 p 值阈值；仅用于初步筛选，后续还会 BH 校正
  pAdjustMethod = "BH",          # ⑥ 多重检验校正方法；发文章常用 BH（FDR）
  verbose       = TRUE         # ⑦ 设为 TRUE 可看到下载进度；国内网慢时常卡这里
)
kegg_result <- kegg_gsea@result