```markdown
# 转录组原始counts转化成TPM的原理以及测试代码

1 先对每条转录本的长度做校正（得到每千碱基的读段数，RPK）；  
2 再把所有转录本的 RPK 求和，作为“百万尺度”的分母；  
3 最后把每个转录本的 RPK 除以该分母并乘 10⁶，得到 TPM；

## pipeline工作流程

step1 运行GTF_Prcess（参考基因组文件处理）得到geneid_efflen  
step2 读取counts文件并进行相应处理  
step3 将counts文件和geneid_efflen文件进行ID匹配  
step4 运行counts2TPM函数  
step5 微调格式并保存

## 3个参考基因组下载网站

https://www.ncbi.nlm.nih.gov/datasets/genome/#!/overview/  
https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/  
https://www.gencodegenes.org/human/release_38.html
```
