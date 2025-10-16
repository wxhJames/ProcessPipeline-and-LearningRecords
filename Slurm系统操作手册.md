
```markdown
# Slurm 超算系统速查手册

> 登录即生产，提交即后台，日志即真相。

---

## 1. 登录
```bash
ssh cn01          # 先跳到计算节点再干活
```

---

## 2. 资源概览
| 命令 | 说明 |
|---|---|
| `sinfo` | 可看队列与节点状态：<br>idle=空，alloc=占，comp=释放中，inact=不可用 |
| `sinfo -n cn01` | 只看 cn01 |
| `sinfo -p cn` | 只看 cn 队列 |

---

## 3. 作业状态
| 命令 | 说明 |
|---|---|
| `squeue` | 我的作业列表 |
| `squeue -j 123456` | 看 123456 号作业 |
| `squeue -u $USER` | 看我自己的全部作业 |
| `squeue -p cn` | 看 cn 队列作业 |

作业状态：R=运行，PD=排队，CG=退出中，S=被挂起。

---

## 4. 交互提交（易断）
```bash
srun -p cn -N 2 -n 4 -t 20 hostname
```
| 选项 | 含义 |
|---|---|
| `-N 2` | 2 节点 |
| `-n 4` | 4 进程 |
| `-c 20` | 每进程 20 核 |
| `-t 20` | 限 20 min |
| `-w cn[01-02]` | 指定节点 |
| `-x cn[01-02]` | 排除节点 |
| `-o out.log -e err.log` | 重定向输出 |

---

## 5. 后台提交（推荐）
1. 写脚本 `job.sh`
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 60
#SBATCH -p cn
srun ./a.out
```
2. 提交
```bash
sbatch job.sh          # 生成 slurm-<JOBID>.out
tail -f slurm-*.out    # 实时看日志
```

---

## 6. 作业管理
| 目的 | 命令 |
|---|---|
| 取消作业 | `scancel <JOBID>` |
| 取消同名 | `scancel -n NAME` |
| 取消队列 | `scancel -p cn` |
| 取消排队 | `scancel -t PENDING` |
| 详情查看 | `scontrol show job <JOBID>` |
| 历史账单 | `sacct -u $USER --format=jobid,partition,state,elapsed` |

---

## 7. 文件 & 编辑
```bash
touch file
vim file
rm file
```
vim 保命三连：  
`Esc` → `:wq` 保存退出；只读时 `:q!` 强制退出。

---

## 8. Conda 环境
```bash
conda create -n torch python=3.10
conda activate torch
```

---

## 9. Python 作业模板
**test.py**
```python
import STCAT, scanpy as sc, os, warnings
warnings.filterwarnings("ignore")

adata = sc.read("data/GSE195832_adata_raw.h5ad")
out   = STCAT.STCAT(adata)

out.write_h5ad("data/GSE195832_adata_raw_processed.h5ad")
out.obs.to_csv("data/GSE195832_adata_raw_processed_obs.csv")
```

**test.sh**
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -p fat
srun python test.py
```
提交：
```bash
sbatch test.sh
tail -f slurm-*.out
```

---

一句话总结：  
**登录→写脚本→sbatch→tail 日志→scancel 不爽就删。**
```