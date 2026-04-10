# run_MrBayes
# 🧬 Bayesian Phylogenetic Tree Pipeline (MrBayes)

本流程用于基于多基因序列构建贝叶斯进化树（Bayesian phylogeny），适用于 DNA 序列数据。

---

## 📦 环境安装

使用 mamba 安装 MrBayes：

```bash
mamba create -n mb -c bioconda mrbayes
```

---

## ⚠️ 注意事项（非常重要）

* 所有序列 ID 严禁包含以下字符：

  * 空格
  * 破折号 `-`
  * 逗号 `,`
  * 点号 `.`
* 必须统一使用 `_` 连接，否则程序会报错

---

## 🚀 分析流程

### 1️⃣ 序列 ID 标准化

```bash
python clean_fasta_ids.py --seqtype dna *.fasta
rm *fasta
```

输出文件：

```
*_clean.fas
```

---

### 2️⃣ 多基因序列合并

```bash
python concat_genes.py -i *_clean.fas -o merged.fasta
```

输出文件：

```
merged.fasta
```

---

### 3️⃣ 多序列比对（MSA）

```bash
mafft --thread 8 --maxiterate 1000 --localpair merged.fasta > merged.aln
```

输出文件：

```
merged.aln
```

---

### 4️⃣ 转换为 MrBayes 输入格式

```bash
python aln_to_mrbayes_nex.py \
-i merged.aln \
-o merged.nex \
--outgroup CP062074_1_Bacillus_velezensis_strain_BSC16a \
--ngen 100000 \
--samplefreq 200
```

参数说明：

| 参数             | 含义                           |
| -------------- | ---------------------------- |
| `--outgroup`   | 外群序列 ID（需存在于 merged.fasta 中） |
| `--ngen`       | MCMC 迭代次数                    |
| `--samplefreq` | 采样频率                         |

输出文件：

```
merged.nex
```

---

### 5️⃣ 运行 MrBayes

```bash
conda activate mb
mpirun -np 8 mb-mpi merged.nex
conda deactivate
```

---

## 📊 输出结果

运行完成后生成：

```
merged.nex.con.tre
```

---

## 🌳 进化树可视化

推荐使用 iTOL 在线工具：

https://itol.embl.de/upload.cgi

上传文件：

```
merged.nex.con.tre
```

即可进行树的可视化与美化。

---

## 📁 文件结构示意

```
project/
│
├── *.fasta
├── *_clean.fas
├── merged.fasta
├── merged.aln
├── merged.nex
├── merged.nex.con.tre
└── scripts/
    ├── clean_fasta_ids.py
    ├── concat_genes.py
    └── aln_to_mrbayes_nex.py
```

---

## 💡 小提示

* 建议 `--ngen` ≥ 100000（复杂数据建议更高）
* 可以在 MrBayes 中检查：

  * ESS
  * 收敛性（PSRF ≈ 1）

---
