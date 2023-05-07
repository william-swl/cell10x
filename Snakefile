import os,time,shutil
pip_dir='/home/william/pipeline/cell10x'
configfile: f'{pip_dir}/sample_config/test.yaml'


# software
cellranger = config['cellranger']

# structure
indir=config['indir']
outdir=config['outdir']
Lsample=config['Lsample']
Dresources=config['Dresources']
Plog=f'{outdir}/0log'
Pstat=f'{outdir}/0stat'
Pqc=f'{outdir}/1qc'
Pcount=f'{outdir}/2count'

for p in [Pqc,Pcount]:
    os.makedirs(p, exist_ok=True)

onstart:
    for r in workflow.rules:
        os.makedirs(f'{Plog}/{r}', exist_ok=True)
    global time_start, smk_file_name
    time_start = time.strftime("%y%m%d_%H%M%S", time.localtime())
    smk_file_name = workflow.snakefile.split('/')[-1]
    shutil.copyfile(workflow.snakefile, f'{Plog}/all/{time_start}_{smk_file_name}')


rule all:
    input:
        # qc
        qc_html = expand(Pqc + '/{sample}/multiqc/{sample}_multiqc_report.html',sample=Lsample),
        # count
        mRNA_count_dir = [Pcount + f'/{sample}/{sample}-mRNA/outs' for sample in Lsample if config[sample]['mRNA']],
        VDJB_count_dir = [Pcount + f'/{sample}/{sample}-VDJB/outs' for sample in Lsample if config[sample]['VDJB']],
        VDJT_count_dir = [Pcount + f'/{sample}/{sample}-VDJT/outs' for sample in Lsample if config[sample]['VDJT']],


rule qc:
    input: indir + '/{sample}',
    output:
        qc_html = Pqc + '/{sample}/multiqc/{sample}_multiqc_report.html',
        qc_data = Pqc + '/{sample}/multiqc/{sample}_multiqc_report_data/multiqc_data.json'
    log: e = Plog + '/qc/{sample}.e', o = Plog + '/qc/{sample}.o'
    benchmark: Plog + '/qc/{sample}.bmk'
    resources: cpus=Dresources['qc_cpus']
    conda: f'{pip_dir}/envs/upstream.yaml'
    shell:"""
        find {input} -name "*.fastq.gz" -or -name "*.fq.gz" | xargs fastqc -o {Pqc}/{wildcards.sample} -t {resources.cpus} 1>{log.o} 2>{log.e}
        multiqc -i {wildcards.sample} -f -o {Pqc}/{wildcards.sample}/multiqc {Pqc}/{wildcards.sample} 1>>{log.o} 2>>{log.e}
        """


##################################
### mRNA
##################################
if lambda wildcards:config[wildcards.sample]['mRNA']:
    rule mRNA_count:
        input: fq_dir = indir + '/{sample}'
        output: mRNA_count_dir = directory(Pcount + '/{sample}/{sample}-mRNA/outs')
        params: gex_ref = config['gex_ref']
        log: e = Plog + '/mRNA_count/{sample}.e', o = Plog + '/mRNA_count/{sample}.o'
        benchmark: Plog + '/mRNA_count/{sample}.bmk'
        resources: cpus=Dresources['mRNA_count_cpus'], count_local_mem=Dresources['count_local_mem']
        conda: f'{pip_dir}/envs/upstream.yaml'
        shell:"""
            cd {Pcount}/{wildcards.sample}
            rm -r {Pcount}/{wildcards.sample}/{wildcards.sample}-mRNA
            {cellranger} count --id={wildcards.sample}-mRNA \\
                --transcriptome={params.gex_ref} \\
                --fastqs={input.fq_dir} \\
                --sample={wildcards.sample}-mRNA \\
                --localcores {resources.cpus} --mempercore 4 \\
                # only use localmem if test on local machine with limited memory 
                --localmem={resources.count_local_mem} \\
                1>{log.o} 2>{log.e}
            cd -
            """

##################################
### VDJB
##################################
if lambda wildcards:config[wildcards.sample]['VDJB']:
    rule VDJB_count:
        input: fq_dir = indir + '/{sample}'
        output: VDJB_count_dir = directory(Pcount + '/{sample}/{sample}-VDJB/outs')
        params: vdj_ref = config['vdj_ref']
        log: e = Plog + '/VDJB_count/{sample}.e', o = Plog + '/VDJB_count/{sample}.o'
        benchmark: Plog + '/VDJB_count/{sample}.bmk'
        resources: cpus=Dresources['VDJB_count_cpus'], count_local_mem=Dresources['count_local_mem']
        conda: f'{pip_dir}/envs/upstream.yaml'
        shell:"""
            cd {Pcount}/{wildcards.sample}
            rm -r {Pcount}/{wildcards.sample}/{wildcards.sample}-VDJB
            {cellranger} vdj --id={wildcards.sample}-VDJB \\
                --reference={params.vdj_ref} \\
                --fastqs={input.fq_dir} \\
                --sample={wildcards.sample}-VDJB \\
                --localcores {resources.cpus} --mempercore 4 \\
                # only use if test on local machine with limited memory 
                --localmem={resources.count_local_mem} \\
                1>{log.o} 2>{log.e}
            cd -
            """


##################################
### VDJT
##################################
if lambda wildcards:config[wildcards.sample]['VDJT']:
    rule VDJT_count:
        input: fq_dir = indir + '/{sample}'
        output: VDJT_count_dir = directory(Pcount + '/{sample}/{sample}-VDJT/outs')
        params: vdj_ref = config['vdj_ref']
        log: e = Plog + '/VDJT_count/{sample}.e', o = Plog + '/VDJT_count/{sample}.o'
        benchmark: Plog + '/VDJT_count/{sample}.bmk'
        resources: cpus=Dresources['VDJT_count_cpus'], count_local_mem=Dresources['count_local_mem']
        conda: f'{pip_dir}/envs/upstream.yaml'
        shell:"""
            cd {Pcount}/{wildcards.sample}
            rm -r {Pcount}/{wildcards.sample}/{wildcards.sample}-VDJT
            {cellranger} vdj --id={wildcards.sample}-VDJT \\
                --reference={params.vdj_ref} \\
                --fastqs={input.fq_dir} \\
                --sample={wildcards.sample}-VDJT \\
                --localcores {resources.cpus} --mempercore 4 \\
                 # only use if test on local machine with limited memory 
                --localmem={resources.count_local_mem} \\
                1>{log.o} 2>{log.e}
            cd -
            """
