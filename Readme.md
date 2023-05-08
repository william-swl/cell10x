# 安装

1. 安装 `conda`和 `mamba`，设置合适的 `channel`
2. 创建运行pipeline的新环境，并安装 `snakemake`

```shell
mamba create -n cell10x
mamba activate cell10x
mamba install snakemake
```

3. 安装pipeline所需依赖

```shell
make setup
```

# 输入

- 将rawdata整理为如下目录结构

```
rawdata/
 |_ test/
   |_ test-mRNA_S...fastq.gz
   |_ test-VDJT_S...fastq.gz
   |_ test-VDJB_S...fastq.gz
 |_ sample1/
   |_ sample1-mRNA...fastq.gz
   |_ sample1-VDJB...fastq.gz
   |_ sample1-FB...fastq.gz

```

# qc

- 对 ` {indir}/{sample}/` 目录下的所有 `fastq.gz `或 `fq.gz`文件产生质控报告，涉及 `{sample}-mRNA_S..., {sample}-VDJB_S...`等

# count


# parse

## mRNA
- 使用SingleR进行细胞类型注释，reference来自[celldex](http://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#34_Monaco_immune_data)
- 默认包含四个注释索引，根据`sample_config`中的`species, celltype_ref_tissue`选定用哪一个
1. 人类通用：human common, HumanPrimaryCellAtlasData
2. 人类免疫：human immune, MonacoImmuneData
3. 小鼠通用：mouse common, MouseRNAseqData
4. 小鼠免疫：mouse immune, ImmGenData
- 输出`Seurat`对象的`rds`，和其`meta.data`的`csv`


## VDJ





# 常见问题

## 软件版本

- `snakemake`：7.25.2

## 安装

`envs/`下的 `xxx.yaml`指定了环境中需要通过 `conda/mamba`安装的程序，与之对应的可能有一个 `xxx.post-deploy.sh`，会在 `conda`环境安装完毕后执行

当出现如下提示时，证明创建的 `conda`环境安装在 `conda_envs/dab33f3628bbe1648ff335ed9bba215f_`，此时 `xxx.post-deploy.sh`正在执行，可通过查看 `conda_envs/dab33f3628bbe1648ff335ed9bba215f_/deploy.log`获取 `xxx.post-deploy.sh`的运行日志

```
Creating conda environment envs/upstream.yaml...
Downloading and installing remote packages.
Running post-deploy script conda_envs/dab33f3628bbe1648ff335ed9bba215f_.post-deploy.sh...
```
