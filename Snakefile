import os,time,shutil,json
pip_dir = os.getcwd()
configfile: f'{pip_dir}/sample_config/test.yaml'


# argument
species = config['species']

# software
cellranger = config['cellranger']

# structure
indir=config['indir']
outdir=config['outdir']
Lsample=config['Lsample']
Plog=f'{outdir}/0log'
Pstat=f'{outdir}/0stat'
Pqc=f'{outdir}/1qc'
Pcount=f'{outdir}/2count'
PmRNA=f'{outdir}/3parse/mRNA'
PFB=f'{outdir}/3parse/FB'
PVDJB=f'{outdir}/3parse/VDJB'
PVDJT=f'{outdir}/3parse/VDJT'
Pfilter=f'{outdir}/4filter'
Pvisualize=f'{outdir}/5visualize'



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
        FB_ref = [Pcount + f'/{sample}/feature_ref.csv' for sample in Lsample if config[sample]['FB']],
        FB_count_dir = [Pcount + f'/{sample}/{sample}-FB/outs' for sample in Lsample if config[sample]['FB']],
        VDJB_count_dir = [Pcount + f'/{sample}/{sample}-VDJB/outs' for sample in Lsample if config[sample]['VDJB']],
        VDJT_count_dir = [Pcount + f'/{sample}/{sample}-VDJT/outs' for sample in Lsample if config[sample]['VDJT']],
        # parse mRNA
        mRNA_csv = [PmRNA + f'/{sample}/mRNA.csv' for sample in Lsample if config[sample]['mRNA']],
        # parse VDJB
        VDJB_csv = [PVDJB + f'/{sample}/VDJB.csv' for sample in Lsample if config[sample]['VDJB']],
        # VDJB_tree = [PVDJB + f'/{sample}/changeo_clone-pass_germ-pass_igphyml-pass.tab' for sample in Lsample if config[sample]['VDJB']],
        # parse VDJT
        VDJT_csv = [PVDJT + f'/{sample}/VDJT.csv' for sample in Lsample if config[sample]['VDJT']],
        # parse FB
        FB_csv = [PFB + f'/{sample}/FB.csv' for sample in Lsample if config[sample]['FB']],
        # filter
        filter_stat = expand(Pfilter + '/{sample}/filter_stat.yaml', sample=Lsample),
        # visualzie
        visualize_html = expand(Pvisualize + '/{sample}.html', sample=Lsample)

rule notebook_init:
    input: 
        mRNA_parse_r='scripts/mRNA_parse.r.ipynb',
        FB_parse_r='scripts/FB_parse.r.ipynb',
        VDJB_parse_r='scripts/VDJB_parse.r.ipynb',
        VDJT_parse_r='scripts/VDJT_parse.r.ipynb',
        filter_r='scripts/filter.r.ipynb',
        visualize_r='scripts/visualize.r.ipynb',
    output: 
        mRNA_parse_r=Plog + '/mRNA_parse.r.ipynb',
        FB_parse_r=Plog + '/FB_parse.r.ipynb',
        VDJB_parse_r=Plog + '/VDJB_parse.r.ipynb',
        VDJT_parse_r=Plog + '/VDJT_parse.r.ipynb',
        filter_r=Plog + '/filter.r.ipynb',
        visualize_r=Plog + '/visualize.r.ipynb',
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
    resources: cpus=config['qc_cpus']
    params: stat_dir = Pstat + '/{sample}/qc'
    conda: f'{pip_dir}/envs/upstream.yaml'
    shell:"""
        find {input} -name "*.fastq.gz" -or -name "*.fq.gz" | xargs fastqc -o {Pqc}/{wildcards.sample} -t {resources.cpus} 1>{log.o} 2>{log.e}
        multiqc -i {wildcards.sample} -f -o {Pqc}/{wildcards.sample}/multiqc {Pqc}/{wildcards.sample} 1>>{log.o} 2>>{log.e}
        mkdir -p {params.stat_dir}
        cp {output.qc_html} {params.stat_dir}
        cp {output.qc_data} {params.stat_dir}
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
        resources: cpus=config['mRNA_count_cpus']
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
            mRNA_stat = PmRNA + '/{sample}/mRNA_stat.yaml',
            stat_dir = directory(Pstat + '/{sample}/mRNA')
        log: notebook = Plog + '/mRNA_parse/{sample}.r.ipynb', e = Plog + '/mRNA_parse/{sample}.e', o = Plog + '/mRNA_parse/{sample}.o'
        benchmark: Plog + '/mRNA_parse/{sample}.bmk'
        resources: cpus=config['mRNA_parse_cpus']
        conda: f'{pip_dir}/envs/RNA.yaml'
        notebook: rules.notebook_init.output.mRNA_parse_r

##################################
### FB
##################################
if lambda wildcards:config[wildcards.sample]['FB']:
    rule FB_ref:
        input: fq_dir = indir + '/{sample}'
        output: 
            FB_ref = Pcount + '/{sample}/feature_ref.csv',
            FB_lib = Pcount + '/{sample}/feature_lib.csv',
        resources: cpus=1
        log: e = Plog + '/FB_ref/{sample}.e', o = Plog + '/FB_ref/{sample}.o'
        run:
            Dsample = config[wildcards.sample]
            with open(output['FB_ref'], 'w') as Ffb:
                Ffb.write('id,name,read,pattern,sequence,feature_type\n')
                for t, d in Dsample['id2seq'].items():
                    if d:
                        for id, seq in d.items():
                            if re.search('[^\x00-\x7F]| ', id):
                                print(f'blank or non-ASCII in feature_barcode id: {id}')
                                raise RuntimeError(f'blank or non-ASCII in feature_barcode id: {id}')
                            Lline = [id,f'{t}_{id}',config['FB_read'],config['FB_pattern'],seq,config['FB_type']]
                            Ffb.write(','.join(Lline) + '\n')
            with open(output['FB_lib'], 'w') as Flib:
                Flib.write('fastqs,sample,library_type\n')
                Lline = [input['fq_dir'], wildcards.sample + '-FB',config['FB_type']]
                Flib.write(','.join(Lline) + '\n')


    rule FB_count:
        input: 
            FB_ref = rules.FB_ref.output.FB_ref,
            FB_lib = rules.FB_ref.output.FB_lib,
        output: FB_count_dir = directory(Pcount + '/{sample}/{sample}-FB/outs')
        log: e = Plog + '/FB_count/{sample}.e', o = Plog + '/FB_count/{sample}.o'
        benchmark: Plog + '/FB_count/{sample}.bmk'
        resources: cpus=config['FB_count_cpus']
        params: gex_ref = config['gex_ref']
        conda: f'{pip_dir}/envs/upstream.yaml'
        shell:"""
            cd {Pcount}/{wildcards.sample}
            rm -r {Pcount}/{wildcards.sample}/{wildcards.sample}-FB
            {cellranger} count --id={wildcards.sample}-FB \\
                --transcriptome={params.gex_ref} \\
                --feature-ref {input.FB_ref} \\
                --libraries {input.FB_lib} \\
                --localcores {resources.cpus} --mempercore 4 \\
                1>{log.o} 2>{log.e}
            cd -
            """

    rule FB_parse:
        input:
            FB_count_dir = rules.FB_count.output.FB_count_dir,
            FB_parse_r = rules.notebook_init.output.FB_parse_r
        output:
            FB_csv = PFB + '/{sample}/FB.csv',
            FB_stat = PFB + '/{sample}/FB_stat.yaml',
            stat_dir = directory(Pstat + '/{sample}/FB')
        log: notebook = Plog + '/FB_parse/{sample}.r.ipynb', e = Plog + '/FB_parse/{sample}.e', o = Plog + '/FB_parse/{sample}.o'
        benchmark: Plog + '/FB_parse/{sample}.bmk'
        resources: cpus=config['FB_parse_cpus']
        conda: f'{pip_dir}/envs/RNA.yaml'
        notebook: rules.notebook_init.output.FB_parse_r


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
        resources: cpus=config['VDJB_count_cpus']
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

    rule VDJB_igblast:
        input: VDJB_count_dir = rules.VDJB_count.output.VDJB_count_dir
        output:
            VDJB_orf_nt_fa = PVDJB + '/{sample}/seq_orf_nt.fasta',
            VDJB_igblast_tsv = PVDJB + '/{sample}/igblast_airr.tsv',
            VDJB_igblast_txt = PVDJB + '/{sample}/igblast_blast.txt',
        log: e = Plog + '/VDJB_igblast/{sample}.e', o = Plog + '/VDJB_igblast/{sample}.o'
        benchmark: Plog + '/VDJB_igblast/{sample}.bmk'
        resources: cpus=config['VDJB_igblast_cpus']
        params: 
            igblast_VDJB_ref_prefix = config['igblast_VDJB_ref_prefix'],
            igblast_aux = config['igblast_aux']
        conda: f'{pip_dir}/envs/VDJ.yaml'
        shell:"""
            # fetch sequences in ORF
            R -e " \
                x <- readr::read_tsv('{input.VDJB_count_dir}/airr_rearrangement.tsv'); \
                genogamesh::parse_CellRanger_vdjseq(x, file='{output.VDJB_orf_nt_fa}', fa_content='seq_orf_nt') \
                " 1>>{log.o} 2>>{log.e}
            
            # airr format
            igblastn -organism {species} -germline_db_V {params.igblast_VDJB_ref_prefix}V -domain_system imgt\\
                -germline_db_D {params.igblast_VDJB_ref_prefix}D -germline_db_J {params.igblast_VDJB_ref_prefix}J \\
                -auxiliary_data {params.igblast_aux} \\
                -show_translation -outfmt 19 -num_threads {resources.cpus} \\
                -query {output.VDJB_orf_nt_fa} \\
                -out {output.VDJB_igblast_tsv} 1>>{log.o} 2>>{log.e}
            
            # blast format
            igblastn -organism {species} -germline_db_V {params.igblast_VDJB_ref_prefix}V -domain_system imgt\\
                -germline_db_D {params.igblast_VDJB_ref_prefix}D -germline_db_J {params.igblast_VDJB_ref_prefix}J \\
                -auxiliary_data {params.igblast_aux} \\
                -show_translation -outfmt '7 std qseq sseq btop' -num_threads {resources.cpus} \\
                -query {output.VDJB_orf_nt_fa} \\
                -out {output.VDJB_igblast_txt} 1>>{log.o} 2>>{log.e}
            """

    rule VDJB_changeo:
        input:
            VDJB_orf_nt_fa = rules.VDJB_igblast.output.VDJB_orf_nt_fa,
            VDJB_count_dir = rules.VDJB_count.output.VDJB_count_dir,
            VDJB_igblast_txt = rules.VDJB_igblast.output.VDJB_igblast_txt
        output: 
            VDJB_changeo_db = PVDJB + '/{sample}/changeo_db-pass.tsv',
            VDJB_changeo = PVDJB + '/{sample}/changeo_clone-pass.tsv',
            VDJB_changeo_fail = PVDJB + '/{sample}/changeo_clone-fail.tsv',
        log: e = Plog + '/VDJB_changeo/{sample}.e', o = Plog + '/VDJB_changeo/{sample}.o'
        benchmark: Plog + '/VDJB_changeo/{sample}.bmk'
        resources: cpus=config['VDJB_changeo_cpus']
        conda: f'{pip_dir}/envs/VDJ.yaml'
        params: 
            changeo_VB_ref = config['changeo_VB_ref'],
            changeo_DB_ref = config['changeo_DB_ref'],
            changeo_JB_ref = config['changeo_JB_ref'],
            outdir = PVDJB + '/{sample}'
        shell:"""
            MakeDb.py igblast -i {input.VDJB_igblast_txt} \\
                -r {params.changeo_VB_ref} {params.changeo_DB_ref} {params.changeo_JB_ref} \\
                -s {input.VDJB_orf_nt_fa} --outname changeo --regions default \\
                --failed --partial --format airr --extended 1>>{log.o} 2>>{log.e}
            
            DefineClones.py -d {output.VDJB_changeo_db} --failed --act set --nproc {resources.cpus}\\
                --outname changeo --model ham --norm len --dist 0.15 1>>{log.o} 2>>{log.e}

            # in case all contigs are failed
            touch {output.VDJB_changeo}
            """

    rule VDJB_tree:
        input:
            VDJB_changeo = rules.VDJB_changeo.output.VDJB_changeo,
        output:
            VDJB_changeo_gm  = PVDJB + '/{sample}/changeo_clone-pass_germ-pass.tsv',
            VDJB_tree = PVDJB + '/{sample}/changeo_clone-pass_germ-pass_igphyml-pass.tab',
        log: e = Plog + '/VDJB_tree/{sample}.e', o = Plog + '/VDJB_tree/{sample}.o'
        benchmark: Plog + '/VDJB_tree/{sample}.bmk'
        resources: cpus=config['VDJB_tree_cpus']
        conda: f'{pip_dir}/envs/VDJ.yaml'
        params: 
            changeo_VB_ref = config['changeo_VB_ref'],
            changeo_DB_ref = config['changeo_DB_ref'],
            changeo_JB_ref = config['changeo_JB_ref'],
            outdir = PVDJB + '/{sample}'
        shell:"""
            CreateGermlines.py -d {input.VDJB_changeo} -g dmask --cloned \\
                -r {params.changeo_VB_ref} {params.changeo_DB_ref} {params.changeo_JB_ref} \\
                1>>{log.o} 2>>{log.e}

            BuildTrees.py -d {output.VDJB_changeo_gm} --collapse \\
                --sample 3000 --igphyml --clean all --nproc {resources.cpus} \\
                1>>{log.o} 2>>{log.e}
        """

    rule VDJB_anarci:
        input: VDJB_count_dir = rules.VDJB_count.output.VDJB_count_dir
        output:
            VDJB_orf_aa_fa = PVDJB + '/{sample}/seq_orf_aa.fasta',
            VDJB_anarci_H = PVDJB + '/{sample}/anarci_H.csv',
            VDJB_anarci_KL = PVDJB + '/{sample}/anarci_KL.csv',
        log: e = Plog + '/VDJB_anarci/{sample}.e', o = Plog + '/VDJB_anarci/{sample}.o'
        benchmark: Plog + '/VDJB_anarci/{sample}.bmk'
        resources: cpus=config['VDJB_anarci_cpus']
        params: outdir = PVDJB + '/{sample}', ext_numbering = config['ext_numbering']
        conda: f'{pip_dir}/envs/VDJ.yaml'
        shell:"""
            # fetch sequences in ORF
            R -e " \
                x <- readr::read_tsv('{input.VDJB_count_dir}/airr_rearrangement.tsv'); \
                genogamesh::parse_CellRanger_vdjseq(x, file='{output.VDJB_orf_aa_fa}', fa_content='seq_orf_aa') \
                " 1>>{log.o} 2>>{log.e}
            
            # airr format
            cd {params.outdir}
            ANARCI -i {output.VDJB_orf_aa_fa} -o anarci -ht anarci_hittable.txt \\
                --use_species {species} --restrict ig -s {params.ext_numbering} --csv --ncpu {resources.cpus} \\
                --assign_germline 1>>{log.o} 2>>{log.e}
            """


    rule VDJB_parse:
        input:
            VDJB_count_dir = rules.VDJB_count.output.VDJB_count_dir,
            VDJB_igblast_tsv = rules.VDJB_igblast.output.VDJB_igblast_tsv,
            VDJB_igblast_txt = rules.VDJB_igblast.output.VDJB_igblast_txt,
            VDJB_changeo = rules.VDJB_changeo.output.VDJB_changeo,
            VDJB_changeo_fail = rules.VDJB_changeo.output.VDJB_changeo_fail,
            VDJB_anarci_H = rules.VDJB_anarci.output.VDJB_anarci_H,
            VDJB_anarci_KL = rules.VDJB_anarci.output.VDJB_anarci_KL,
            VDJB_parse_r = rules.notebook_init.output.VDJB_parse_r
        output:
            VDJB_csv = PVDJB + '/{sample}/VDJB.csv',
            VDJB_stat = PVDJB + '/{sample}/VDJB_stat.yaml',
            stat_dir = directory(Pstat + '/{sample}/VDJB')
        log: notebook = Plog + '/VDJB_parse/{sample}.r.ipynb', e = Plog + '/VDJB_parse/{sample}.e', o = Plog + '/VDJB_parse/{sample}.o'
        benchmark: Plog + '/VDJB_parse/{sample}.bmk'
        resources: cpus=config['VDJB_parse_cpus']
        conda: f'{pip_dir}/envs/VDJ.yaml'
        notebook: rules.notebook_init.output.VDJB_parse_r


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
        resources: cpus=config['VDJT_count_cpus']
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

    rule VDJT_parse:
        input:
            VDJT_count_dir = rules.VDJT_count.output.VDJT_count_dir,
            VDJT_parse_r = rules.notebook_init.output.VDJT_parse_r
        output:
            VDJT_csv = PVDJT + '/{sample}/VDJT.csv',
            VDJT_stat = PVDJT + '/{sample}/VDJT_stat.yaml',
            stat_dir = directory(Pstat + '/{sample}/VDJT')
        log: notebook = Plog + '/VDJT_parse/{sample}.r.ipynb', e = Plog + '/VDJT_parse/{sample}.e', o = Plog + '/VDJT_parse/{sample}.o'
        benchmark: Plog + '/VDJT_parse/{sample}.bmk'
        resources: cpus=config['VDJT_parse_cpus']
        conda: f'{pip_dir}/envs/VDJ.yaml'
        notebook: rules.notebook_init.output.VDJT_parse_r

##################################
### filter
##################################
def get_filter_input(wildcards):
    dirs = []
    if config[wildcards.sample]['mRNA']:
        dirs.append(rules.mRNA_parse.output.stat_dir)
    if config[wildcards.sample]['FB']:
        dirs.append(rules.FB_parse.output.stat_dir)
    if config[wildcards.sample]['VDJB']:
        dirs.append(rules.VDJB_parse.output.stat_dir)
    if config[wildcards.sample]['VDJT']:
        dirs.append(rules.VDJT_parse.output.stat_dir)
    return dirs

rule filter:
    input:
        filter_input = get_filter_input,
        filter_r = rules.notebook_init.output.filter_r
    output:
        filter_dir = directory(Pfilter + '/{sample}'),
        filter_stat = Pfilter + '/{sample}/filter_stat.yaml'
    params: stat_dir = Pstat + '/{sample}'
    log: notebook = Plog + '/filter/{sample}.r.ipynb', e = Plog + '/filter/{sample}.e', o = Plog + '/filter/{sample}.o'
    benchmark: Plog + '/filter/{sample}.bmk'
    resources: cpus=config['filter_cpus']
    conda: f'{pip_dir}/envs/visualize.yaml'
    notebook: rules.notebook_init.output.filter_r

##################################
### visualize
##################################

rule visualize:
    input:
        filter_dir = rules.filter.output.filter_dir
    output:
        visualize_rds = Pvisualize + '/{sample}.rds'
    log: notebook = Plog + '/visualize/{sample}.r.ipynb', e = Plog + '/visualize/{sample}.e', o = Plog + '/visualize/{sample}.o'
    params: 
        stat_dir = Pstat + '/{sample}',
        echarts_theme = f'{pip_dir}/src/echarts_theme/mytheme.json',
        metadata = config['metadata']
    benchmark: Plog + '/visualize/{sample}.bmk'
    resources: cpus=config['visualize_cpus']
    conda: f'{pip_dir}/envs/visualize.yaml'
    notebook: rules.notebook_init.output.visualize_r

rule visualize_rmd:
    input:
        visualize_rds = rules.visualize.output.visualize_rds
    output:
        visualize_html = Pvisualize + '/{sample}.html'
    log: e = Plog + '/visualize_rmd/{sample}.e', o = Plog + '/visualize_rmd/{sample}.o'
    params: css_dir = f'{pip_dir}/src/visualize_css'
    benchmark: Plog + '/visualize_rmd/{sample}.bmk'
    resources: cpus=config['visualize_rmd_cpus']
    conda: f'{pip_dir}/envs/visualize.yaml'
    script: f'{pip_dir}/scripts/visualize.Rmd'
