```markdown
# GSEA分析的3个主要函数

1.GSEA()
2.gseKEGG()  经常连不上网
3.gseGO()

# pipeline工作流程

step1 读取数据文件
step2 构建genelist
step3 GSEA分析

# 数据文件格式

经过差异基因分析（limma/deseq）表达得到的数据，一个list，数据格式如图
![图](data_file.png)
