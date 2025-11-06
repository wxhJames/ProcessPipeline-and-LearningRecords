# 数据读取
data_file<-readRDS("tem_file.rds")

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

# 三种函数进行GSEA分析  
# GSEA()
# gseKEGG()  解释“信号通路”
# gseGo()  解释“生物学过程/功能”

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

############################################################################
# h.all.v2025.1.Hs.symbols.gmt  GSEA官网下载

# 读取基因集
# 列名：TERM, GENE（都是 SYMBOL）
gmt_sym <- read.gmt("h.all.v2025.1.Hs.symbols.gmt") 

# 构建genelist
geneList_Symbol <- data_file %>%
  select(Gene, logFC) %>%      # 只要基因名 + 变化倍数
  deframe() %>%                # 转 named numeric
  sort(decreasing = TRUE)      # 保证降序（GSEA 要求）


## 通用 GSEA 富集分析 ----------------------------------------------------
gsea_res <- GSEA(
  geneList      = geneList_Symbol ,      # ① 已按 log2FC 降序排好的 named numeric 向量；名字须与 gmt_df$GENE 同类型
  TERM2GENE     = gmt_sym,        # ② 两列 data.frame：第 1 列基因集名，第 2 列基因 ID；位置决定含义，列名任意
  minGSSize     = 1,            # ③ 基因集本身 ≥15 个基因才纳入检验，太小统计力不足
  maxGSSize     = 500,           # ④ 基因集本身 >500 个基因就跳过，太泛容易稀释信号
  pvalueCutoff  = 0.05,          # ⑤ 原始 p 值初筛阈值；最终显著性看 p.adjust 列
  pAdjustMethod = "BH",          # ⑥ 多重检验校正方法；BH 即 Benjamini-Hochberg（FDR），发文通用
  nPerm         = 1000,          # ⑦ 置换次数；1000 次是常见精度，越大越准但越慢
  verbose       = TRUE           # ⑧ 打印进度条；国内网慢或批跑时可设 FALSE 减少日志
)

## 返回对象
# gsea_res@result  数据框：每条基因集的 ES、NES、pvalue、p.adjust、geneID（leading edge）
# 可直接喂给 dotplot()/gseaplot2()/enrichmap() 画图

gsea_result <- gsea_res@result