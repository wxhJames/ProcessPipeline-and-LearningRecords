
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
## GO 富集分析（gseGO）-------------------------------------------------
go_gsea <- gseGO(
  geneList      = geneList_ENTREZ,      # ① 已按 log2FC 降序排好的 named numeric 向量；名字必须与 OrgDb 的 keyType 一致（默认 EntrezID）
  OrgDb         = org.Hs.eg.db,  # ② 人类注释数据库；小鼠改用 org.Mm.eg.db，大鼠 org.Rn.eg.db
  ont           = "ALL",         # ③ 同时扫描 BP（生物过程）、CC（细胞组分）、MF（分子功能）；可单独写 "BP" / "CC" / "MF"
  minGSSize     = 10,            # ④ GO 术语本身至少包含 10 个基因，低于此数不参与检验（太小统计力不足）
  maxGSSize     = 500,           # ⑤ GO 术语本身最多 500 个基因，过大容易泛化、稀释信号
  pvalueCutoff  = 0.05,          # ⑥ 原始 p 值筛选阈值；仅做初步过滤，最终看 p.adjust 列
  pAdjustMethod = "BH",          # ⑦ 多重检验校正方法；BH 即 Benjamini-Hochberg（FDR），发文章最常用
  verbose       = TRUE          # ⑧ 是否打印下载/运行日志；国内网慢时可设 TRUE 看进度
)
go_result   <- go_gsea@result
