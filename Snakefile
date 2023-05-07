import os,time,shutil,json
pip_dir = os.getcwd()
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
PmRNA=f'{outdir}/3parse/mRNA'

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
        # parse
        mRNA_csv = [PmRNA + f'/{sample}/mRNA.csv' for sample in Lsample if config[sample]['mRNA']],

rule notebook_init:
    input: 
        mRNA_parse_r='scripts/mRNA_parse.r.ipynb'
    output: 
        mRNA_parse_r=Plog + '/mRNA_parse.r.ipynb'
    resources: cpus=1
    log: e = Plog + '/notebook_init.e', o = Plog + '/notebook_init.o'
    run:
        for k, f in input.items():
            notebook = json.loads(open(f, 'r').read())
            if f.endswith('.r.ipynb'):
                notebook['metadata']['kernelspec']['name'] = 'ir'
                notebook['metadata']['kernelspec']['display_name'] = 'R'
            elif f.endswith('.py.ipynb'):
                notebook['metadata']['kernelspec']['name'] = 'python3'
                notebook['metadata']['kernelspec']['display_name'] = 'Python 3 (ipykernel)'
            with open(output[k], 'w') as out:
                json.dump(notebook, out)


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
        resources: cpus=Dresources['mRNA_count_cpus']
        conda: f'{pip_dir}/envs/upstream.yaml'
        shell:"""
            cd {Pcount}/{wildcards.sample}
            rm -r {Pcount}/{wildcards.sample}/{wildcards.sample}-mRNA
            {cellranger} count --id={wildcards.sample}-mRNA \\
                --transcriptome={params.gex_ref} \\
                --fastqs={input.fq_dir} \\
                --sample={wildcards.sample}-mRNA \\
                --localcores {resources.cpus} --mempercore 4 \\
                1>{log.o} 2>{log.e}
            cd -
            """

    rule mRNA_parse:
        input:
            mRNA_count_dir = rules.mRNA_count.output.mRNA_count_dir,
            mRNA_parse_r = rules.notebook_init.output.mRNA_parse_r
        output:
            mRNA_csv = PmRNA + '/{sample}/mRNA.csv',
            mRNA_rds = PmRNA + '/{sample}/mRNA.rds',
            mRNA_stat = PmRNA + '/{sample}/mRNA_stat.yaml'
        log: notebook = Plog + '/mRNA_parse/{sample}.r.ipynb', e = Plog + '/mRNA_parse/{sample}.e', o = Plog + '/mRNA_parse/{sample}.o'
        benchmark: Plog + '/mRNA_parse/{sample}.bmk'
        resources: cpus=Dresources['mRNA_parse_cpus']
        conda: f'{pip_dir}/envs/RNA.yaml'
        notebook: rules.notebook_init.output.mRNA_parse_r


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
        resources: cpus=Dresources['VDJB_count_cpus']
        conda: f'{pip_dir}/envs/upstream.yaml'
        shell:"""
            cd {Pcount}/{wildcards.sample}
            rm -r {Pcount}/{wildcards.sample}/{wildcards.sample}-VDJB
            {cellranger} vdj --id={wildcards.sample}-VDJB \\
                --reference={params.vdj_ref} \\
                --fastqs={input.fq_dir} \\
                --sample={wildcards.sample}-VDJB \\
                --localcores {resources.cpus} --mempercore 4 \\
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
        resources: cpus=Dresources['VDJT_count_cpus']
        conda: f'{pip_dir}/envs/upstream.yaml'
        shell:"""
            cd {Pcount}/{wildcards.sample}
            rm -r {Pcount}/{wildcards.sample}/{wildcards.sample}-VDJT
            {cellranger} vdj --id={wildcards.sample}-VDJT \\
                --reference={params.vdj_ref} \\
                --fastqs={input.fq_dir} \\
                --sample={wildcards.sample}-VDJT \\
                --localcores {resources.cpus} --mempercore 4 \\
                1>{log.o} 2>{log.e}
            cd -
            """
