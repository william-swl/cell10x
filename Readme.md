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
   |_ test-mRNA
   |_ test-VDJT
   |_ test-VDJB
 |_ sample1/
   |_ sample1-mRNA
   |_ sample1-VDJB
   |_ sample1-FB

```

# qc

- 对 ` {indir}/{sample}/` 目录下的所有 `fastq.gz `或 `fq.gz`文件产生质控报告，涉及 `{sample}-mRNA_S..., {sample}-VDJB_S...`等

# count

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
