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

# 配置

- 在 `sample_config`目录下添加配置文件，文件名为 `{sample}.yaml`，如 `test.yaml`，`batch1.yaml`, `batch2.yaml`等
- 参照 `sample_config/test.yaml`的格式，在 `Key Variables`部分指定该batch的输入输出路径、包含的样本、物种等关键信息
- 在 `Sample Variables`部分，指定每个样本包含的测序文库，如 `mRNA, VDJB, VDJT, FB`等。也可以在此处单独指定每个样本的运行参数，如使用特殊的细胞类型注释索引、特殊的过滤条件等，这会覆盖 `Default Variables`部分与其名称相同的参数
- 在 `Default Variables`部分，指定默认的运行参数
- 在 `Resources`部分，指定运行时占用的cpu资源。每个cpu对应的内存数在 `~/.config/snakemake/cell10x/config.yaml`中指定
- 在 `Constant`部分，指定一般不会改变的常量，例如索引文件的位置

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
1. 不再运行igblast、changeo、ANARCI，仅从`CellRanger`的结果取
2. `seq_nt`不再从起始密码子开始，而是序列首端
3. 由于`CellRanger`的输出就是这样，`seq_nt`和`seq_align_nt_H`完全一致


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
