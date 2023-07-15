# 安装

1. 安装 `conda`和 `mamba`，设置合适的 `channel`
2. 创建运行pipeline的新环境，并安装 `snakemake`

```shell
mamba create -n cell10x
mamba activate cell10x
mamba install snakemake
pip install pycones
```

3. 访问[10x官方](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/7.1/)，获取 `Cellranger`的下载链接

![img](doc/fig/cellranger_download.png)

4. 运行脚本搭建运行环境

```
bash init.sh
```

- 以上命令会：
  - 将 `smk_profiles`内的匿名 `conda`路径更改为当前目录
  - 在 `~./config/snakemake`下创建 `profile`，以供 `snakemake --profile`调用
  - 将 `sample_config/test.yaml`内的测试数据集路径更改为当前目录
  - 安装 `Cellranger`、所需参考数据集
  - 设置 `conda`的channel优先选择 `conda-forge`，以[尽量避免兼容性问题](https://conda-forge.org/docs/user/tipsandtricks.html)
  - 运行 `snakemake --conda-create-envs-only`，创建所需的匿名 `conda`环境
- 如果一切顺利，将不需要额外下载任何资源

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

# 配置

- 在 `sample_config`目录下添加配置文件，文件名为 `{sample}.yaml`，如 `test.yaml`，`batch1.yaml`, `batch2.yaml`等
- 参照 `sample_config/test.yaml`的格式，在 `Key Variables`部分指定该batch的输入输出路径、包含的样本、物种等关键信息
- 在 `Sample Variables`部分，指定每个样本包含的测序文库，如 `mRNA, VDJB, VDJT, FB`等。也可以在此处单独指定每个样本的运行参数，如使用特殊的细胞类型注释索引、特殊的过滤条件等，这会覆盖 `Default Variables`部分与其名称相同的参数
- 在 `Default Variables`部分，指定默认的运行参数
- 在 `Resources`部分，指定运行时占用的cpu资源。每个cpu对应的内存数在 `~/.config/snakemake/cell10x/config.yaml`中指定
- 在 `Constant`部分，指定一般不会改变的常量，例如索引文件的位置

## Feature barcode

- `id2seq`，列出FB名和序列的对应关系
  - `HT：Hashtag`，用于区分细胞身份的标签，例如不同的实验组、不同的样本来源。理想情况下，每个细胞都有且仅有一种HT
  - `NC：Negative control`，阴性对照，例如OVA。理想情况下，细胞上不应该有NC
  - `BD：Binding`，抗原结合
- `FB_corr`，FB矫正系数。根据`id2seq`下FB的名称`fb`通过正则表达式`fb$`寻找相应的列，然后将这些列的FB umi数值乘矫正系数
- `Nlim_HT: 10`。HT如果低于此数值，则标记该细胞为“low hashtag”
- `Nratio_HT: 0.9`。主要HT占全部HT的最低比例，如果低于此数值，则该细胞标记为“mixed”
- `Nratio_BD: 0.2`。某种BD占全部BD的比例，如果超过该数值，则在BD_type中添加此种BD
- `Nratio_NC: 0.1`。NC占BD+NC的比例，如果高于此数值，则会被标记为“NC_load: TRUE”

## filter

- `flt_mode`，过滤模式，形式是 `name: [lib1, lib2, ...]`。先将指定库的未过滤细胞取交集，然后应用所有库的过滤条件。如果mRNA库既有B细胞、也有T细胞，可以设置两个过滤模式分别处理
- `mRNA_gene_flt: 10`：滤除mRNA库基因数低于此数值的细胞
- `mRNA_umi_flt: 10`：滤除mRNA库umi数低于此数值的细胞
- `mRNA_mt_percent_flt: 10`：滤除mRNA库线粒体基因umi数占比高于此数值（百分数）的细胞
- `VDJB_umi_flt: 5`：滤除VDJT库umi数低于此数值的细胞
- `VDJT_umi_flt: 5`：滤除VDJB库umi数低于此数值的细胞
- `FB_umi_flt: 10`：滤除FB库BD umi数低于此数值的细胞

# qc

- 对 ` {indir}/{sample}/` 目录下的所有 `fastq.gz `或 `fq.gz`文件产生质控报告，涉及 `{sample}-mRNA_S..., {sample}-VDJB_S...`等

```shell
1qc
└── test
    ├── multiqc
    │ ├── test_multiqc_report.html
    │ └── test_multiqc_report_data
    │     ├── multiqc.log
    │     ├── multiqc_citations.txt
    │     ├── multiqc_data.json
    │     ├── multiqc_fastqc.txt
    │     ├── multiqc_general_stats.txt
    │     └── multiqc_sources.txt
    ├── test-VDJB_S85_L001_R1_001_fastqc.html
    ├── test-VDJB_S85_L001_R1_001_fastqc.zip
    ...
```

# count

- 使用 `Cellranger`比对、计数

```
2count
└── test
    ├── __test-mRNA.mro
    ├── test-VDJB
    ├── test-VDJT
    └── test-mRNA
        ├── SC_RNA_COUNTER_CS
        │ ├── CELLRANGER_PREFLIGHT
        │ ├── CELLRANGER_PREFLIGHT_LOCAL
        │ ├── FULL_COUNT_INPUTS
        │ ├── GET_AGGREGATE_BARCODES_OUT
        │ ├── SC_MULTI_CORE
        │ ├── WRITE_GENE_INDEX
        │ ├── _STRUCTIFY
        │ └── fork0
        ├── _cmdline
        ├── _filelist
        ├── _finalstate
        ├── _invocation
        ├── _jobmode
        ├── _log
        ├── _mrosource
        ├── _perf
        ├── _sitecheck
        ├── _tags
        ├── _timestamp
        ├── _uuid
        ├── _vdrkill
        ├── _versions
        ├── outs
        │ ├── analysis
        │ ├── cloupe.cloupe
        │ ├── filtered_feature_bc_matrix
        │ ├── filtered_feature_bc_matrix.h5
        │ ├── metrics_summary.csv
        │ ├── molecule_info.h5
        │ ├── possorted_genome_bam.bam
        │ ├── possorted_genome_bam.bam.bai
        │ ├── raw_feature_bc_matrix
        │ ├── raw_feature_bc_matrix.h5
        │ └── web_summary.html
        └── test-mRNA.mri.tgz

```

# parse

## mRNA

- 使用SingleR进行细胞类型注释，reference来自[celldex](http://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#34_Monaco_immune_data)
- 默认包含四个注释索引，根据 `sample_config`中的 `species, celltype_ref_tissue`选定用哪一个

1. 人类通用：human common, HumanPrimaryCellAtlasData
2. 人类免疫：human immune, MonacoImmuneData
3. 小鼠通用：mouse common, MouseRNAseqData
4. 小鼠免疫：mouse immune, ImmGenData

- 输出 `Seurat`对象的 `rds`，和其 `meta.data`的 `csv`

![](doc/fig/mRNA_parse_csv.png)

## VDJ

1. `CellRanger vdj`的输出中，`airr_rearrangement.tsv`可以提取出从起始密码子开始的核酸序列，因此先产生 `seq_orf_nt.fasta`，后续分析基于此序列。最终的 `seq_nt`也取自该文件。 `all_contig_annotations.csv `里有reads和umis信息，留待整合。其他大多注释 `IgBlast`也会给且更详细，最终都采用 `IgBlast`的结果
2. 对 `seq_orf_nt.fasta`使用 `IgBlast`，指定IMGT编号系统，输出 `airr`格式用于获取主要的VDJ注释，输出 `blast`格式用于获取mismatch计算SHM
3. 使用 `Change-O`，输出 `changeo_clone-pass.tsv`用于获取克隆信息
4. 使用 `ANARCI`，输入氨基酸序列，获取额外编号系统下的CDR氨基酸序列（默认为Chothia）。而且由于 `ANARCI`输出的V-domain氨基酸序列比 `IgBlast`输出的更完整，用于结果中的 `seq_align_aa`
5. 使用vdj的总突变数量，除以 `seq_align_nt`的长度，得到shm
6. 对于同一个cell、同一种contig，如果有多条，取uims最多的那条，unique标记为FALSE
7. 根据轻重链分别展开以上信息

VDJT与VDJB的区别：

1. 不再运行igblast、changeo、ANARCI，仅从 `CellRanger`的结果取
2. `seq_nt`不再从起始密码子开始，而是序列首端
3. 由于 `CellRanger`的输出就是这样，`seq_nt`和 `seq_align_nt_H`完全一致

| 项目             | 1                          | 2                      | 3                       | 4                              | 5                      | 6            |
| ---------------- | -------------------------- | ---------------------- | ----------------------- | ------------------------------ | ---------------------- | ------------ |
| 软件             | Cellranger                 | Cellranger             | IgBlast                 | IgBlast                        | Change-O               | ANARCI       |
| 输出文件名       | all_contig_annotations.csv | airr_rearrangement.tsv | -outfmt 19              | -outfmt '7 std qseq sseq btop' | changeo_clone-pass.tsv | anarci_H.csv |
| contig id        | contig_id                  | sequence_id            | sequence_id             | Query                          | sequence_id            | id           |
| contig类型       | chain                      |                        | locus                   |                                | locus                  | chain_type   |
| vdj基因          | v_gene                     | x_call                 |                         |                                |                        |              |
| 详细vdj基因      |                            |                        | x_call                  | √                             | x_call                 | x_gene       |
| c基因            | c_gene                     | c_call                 |                         |                                | c_call                 |              |
| vdj-nt序列       |                            |                        | x_sequence_alignment    |                                |                        |              |
| vdj-aa序列       |                            |                        | x_sequence_alignment_aa |                                |                        | numbering    |
| cdr-aa序列       | cdrx                       |                        | cdrx_aa                 |                                |                        |              |
| cdr-nt序列       | cdrx_nt                    |                        | cdrx                    |                                | cdrx                   |              |
| fwr-aa序列       | fwrx                       |                        | fwrx_aa                 |                                |                        |              |
| fwr-nt序列       | fwrx_nt                    |                        | fwrx                    |                                | fwrx                   |              |
| np序列           |                            |                        | np1/2                   |                                |                        |              |
| reads            | reads                      |                        |                         |                                |                        |              |
| umis             | umis                       |                        |                         |                                |                        |              |
| clone            | raw_clonotype_id           | clone_id               |                         |                                | clone_id               |              |
| 全长aa序列       |                            | sequence_aa            |                         |                                |                        |              |
| 全长nt序列       |                            | sequence               | sequence                |                                | sequence               |              |
| 比对aa序列       |                            |                        | sequence_alignment_aa   |                                |                        |              |
| 比对nt序列       |                            | sequence_alignment     | sequence_alignment      |                                | sequence_alignment     |              |
| 比对germline序列 |                            | germline_alignment     | germline_alignment      |                                | germline_alignment     |              |
| vdj起止位点      |                            | x_sequence_start/end   | x_sequence_start/end    | √                             | x_sequence_start       |              |
| c起止位点        |                            | c_sequence_start       |                         |                                |                        |              |
| cdr起止位点      |                            |                        | cdrx_start/end          |                                |                        |              |
| 错配数           |                            |                        |                         | √                             |                        |              |
| gap数            |                            |                        |                         | √                             |                        |              |

## tree and graph

- 根据`ChangeO`掩蔽d基因后的结果，创建克隆树。每个`ChangeO`克隆将会构建为一棵树，仅包括重链或仅包括轻链
- 根据`igblast`结果中的`sequence_alignment`和`germline_alignment`（被重命名为`seq_align_nt_H`，`gm_align_nt_H`），先用`shazam::collapseClones`生成克隆内公共序列，然后计算不同克隆之间的字符编辑距离
- 使用`MDS`降维，散点图绘制二维关系



# filter

- 用于filter的列先进行初始化，如果值为NA修改为合适的值，使得未通过filter、没有测该库能够区分
- 按照设置的`flt_mode`进行过滤，每种模式都有输出结果
- 对于`Bcell`模式，除了输出过滤后的表格外，也同changeo的输出取一个子集，方便之后B细胞建树
- 将`metadata.csv`中以`#`开头的列添加到输出的表格中

## 可视化

- `visualize`调用 `jupyter notebook`绘制所需图形，并保存到 `.rds`文件中
- `visualize_rmd`使用 `Rmarkdown`，将 `.rds`文件中的图形绘制到 `.html`文件中
- 将`metadata.csv`中不以`#`开头的列，输出到metadata部分
- 存在通过filter，但`clone_changeo_H`未被分配的情况，可视化clone之前滤除了这种情况

[主题编辑器 - Apache ECharts](https://echarts.apache.org/zh/theme-builder.html)

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
