
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