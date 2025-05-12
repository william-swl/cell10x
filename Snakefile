import os,time,shutil,pycones
import pandas as pd
configfile: 'configs/test.yaml'

# argument
species = config['species']

# software
cellranger = os.path.abspath('src/cellranger-7.1.0/bin/cellranger')

# structure
indir=os.path.abspath(config['indir'])
outdir=os.path.abspath(config['outdir'])
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
Ptree=f'{outdir}/6tree'



for p in [Pqc,Pcount]:
    os.makedirs(p, exist_ok=True)

onstart:
    for r in workflow.rules:
        os.makedirs(f'{Plog}/{r}', exist_ok=True)
    global time_start, smk_file_name
    time_start = time.strftime("%y%m%d_%H%M%S", time.localtime())
    smk_file_name = workflow.snakefile.split('/')[-1]
    shutil.copyfile(workflow.snakefile, f'{Plog}/all/{time_start}_{smk_file_name}')
    for f in os.listdir('scripts'):
        if f.endswith('.r.ipynb'):
            pycones.nb_kernel_switch(f'scripts/{f}', f'{Plog}/{f}', kernel='r')
        elif f.endswith('.py.ipynb'):
            pycones.nb_kernel_switch(f'scripts/{f}', f'{Plog}/{f}', kernel='python')


rule all:
    input:
        # qc
        qc_html = expand(Pqc + '/{sample}/multiqc/{sample}_multiqc_report.html',sample=Lsample),
        # count config
        cellranger_config = expand(Pcount + '/{sample}/cellranger_config.csv', sample=Lsample),
        # count
        count_dir = expand(Pcount + '/{sample}/{sample}/outs', sample=Lsample),
        # VDJ
        VDJB_igblast_tsv = [PVDJB + f'/{sample}/igblast_airr.tsv' for sample in Lsample if config[sample]['VDJB']],
        # parse
        mRNA_csv = [PmRNA + f'/{sample}/mRNA.csv' for sample in Lsample if config[sample]['mRNA']],
        VDJB_csv = [PVDJB + f'/{sample}/VDJB.csv' for sample in Lsample if config[sample]['VDJB']],
        VDJT_csv = [PVDJT + f'/{sample}/VDJT.csv' for sample in Lsample if config[sample]['VDJT']],
        FB_csv = [PFB + f'/{sample}/FB.csv' for sample in Lsample if config[sample]['FB']],
        # filter
        filter_stat = expand(Pfilter + '/{sample}/filter_stat.yaml', sample=Lsample),
        # concat stat
        stats = expand(Pstat + '/{sample}/stats.csv', sample=Lsample),
        # batch stat
        batch_stat = Pstat + '/batch_stat.csv' if config['run_batch_stat'] else [],
        # visualzie
        visualize_html = expand(Pvisualize + '/{sample}.html', sample=Lsample),
        # B cell tree
        Bcell_tree_H = [Ptree + f'/{sample}/Bcell_tree_H_igphyml-pass.tab' for sample in Lsample if (config[sample]['VDJB'] and config['run_tree'])]

        

wildcard_constraints:
    sample='[^/]+'


rule qc:
    input: indir + '/{sample}',
    output:
        qc_html = Pqc + '/{sample}/multiqc/{sample}_multiqc_report.html',
        qc_data = Pqc + '/{sample}/multiqc/{sample}_multiqc_report_data/multiqc_data.json'
    log: e = Plog + '/qc/{sample}.e', o = Plog + '/qc/{sample}.o'
    benchmark: Plog + '/qc/{sample}.bmk'
    resources: cpus=config['qc_cpus']
    params: stat_dir = Pstat + '/{sample}/qc'
    conda: 'envs/upstream.yaml'
    shell:"""
        find {input} -name "*.fastq.gz" -or -name "*.fq.gz" | xargs fastqc -o {Pqc}/{wildcards.sample} -t {resources.cpus} 1>{log.o} 2>{log.e}
        multiqc -i {wildcards.sample} -f -o {Pqc}/{wildcards.sample}/multiqc {Pqc}/{wildcards.sample} 1>>{log.o} 2>>{log.e}
        mkdir -p {params.stat_dir}
        cp {output.qc_html} {params.stat_dir}
        cp {output.qc_data} {params.stat_dir}
        """

rule cellranger_config:
    output: cellranger_config = Pcount + '/{sample}/cellranger_config.csv',
    log: notebook = Plog + '/cellranger_config/{sample}.r.ipynb', 
         e = Plog + '/cellranger_config/{sample}.e', o = Plog + '/cellranger_config/{sample}.o'
    benchmark: Plog + '/cellranger_config/{sample}.bmk'
    resources: cpus=1
    params: 
        feature_ref=Pcount + '/{sample}/feature_ref.csv'
    conda: 'envs/RNA.yaml'
    notebook: Plog + '/cellranger_config.r.ipynb'


rule count:
    input: cellranger_config = rules.cellranger_config.output.cellranger_config,
    output: count_dir = directory(Pcount + '/{sample}/{sample}/outs')
    log: e = Plog + '/count/{sample}.e', o = Plog + '/count/{sample}.o'
    benchmark: Plog + '/count/{sample}.bmk'
    resources: cpus=config['count_cpus']
    shell:"""
        rm -r {Pcount}/{wildcards.sample}/{wildcards.sample}
        cd {Pcount}/{wildcards.sample}
        {cellranger} multi --id={wildcards.sample} \\
            --csv={input} --localcores {resources.cpus} --mempercore 4 \\
            1>{log.o} 2>{log.e}
        cd -
        """

##################################
### mRNA
##################################
if lambda wildcards:config[wildcards.sample]['mRNA']:
    rule mRNA_parse:
        input: count_dir = rules.count.output.count_dir
        output:
            mRNA_csv = PmRNA + '/{sample}/mRNA.csv',
            mRNA_rds = PmRNA + '/{sample}/mRNA.rds',
            mRNA_stat = PmRNA + '/{sample}/mRNA_stat.yaml',
            stat_dir = directory(Pstat + '/{sample}/mRNA')
        log: notebook = Plog + '/mRNA_parse/{sample}.r.ipynb', 
             e = Plog + '/mRNA_parse/{sample}.e', o = Plog + '/mRNA_parse/{sample}.o'
        benchmark: Plog + '/mRNA_parse/{sample}.bmk'
        resources: cpus=config['mRNA_parse_cpus']
        conda: 'envs/RNA.yaml'
        notebook: Plog + '/mRNA_parse.r.ipynb'



##################################
### FB
##################################
if lambda wildcards:config[wildcards.sample]['FB']:
    rule FB_parse:
        input: count_dir = rules.count.output.count_dir
        output:
            FB_csv = PFB + '/{sample}/FB.csv',
            FB_stat = PFB + '/{sample}/FB_stat.yaml',
            stat_dir = directory(Pstat + '/{sample}/FB')
        log: notebook = Plog + '/FB_parse/{sample}.r.ipynb', 
             e = Plog + '/FB_parse/{sample}.e', o = Plog + '/FB_parse/{sample}.o'
        benchmark: Plog + '/FB_parse/{sample}.bmk'
        resources: cpus=config['FB_parse_cpus']
        conda: 'envs/RNA.yaml'
        notebook: Plog + '/FB_parse.r.ipynb'


##################################
### VDJB
##################################
if lambda wildcards:config[wildcards.sample]['VDJB']:
    rule VDJB_igblast:
        input: count_dir = rules.count.output.count_dir
        output:
            VDJB_orf_nt_fa = PVDJB + '/{sample}/seq_orf_nt.fasta',
            VDJB_igblast_tsv = PVDJB + '/{sample}/igblast_airr.tsv',
            VDJB_igblast_txt = PVDJB + '/{sample}/igblast_blast.txt',
        log: e = Plog + '/VDJB_igblast/{sample}.e', o = Plog + '/VDJB_igblast/{sample}.o'
        benchmark: Plog + '/VDJB_igblast/{sample}.bmk'
        resources: cpus=config['VDJB_igblast_cpus']
        params: 
            igblast_VDJB_ref_prefix = f'src/igblast-ref/{species}-VDJB/{species}_gl_',
            igblast_aux = f'src/igblast-aux/{species}_gl.aux',
            VDJB_airr = rules.count.output.count_dir + f"/{config['count_VDJB_airr']}"
        conda: 'envs/VDJ.yaml'
        shell:"""
            # fetch sequences in ORF
            R -e " \
                x <- readr::read_tsv('{params.VDJB_airr}'); \
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
            VDJB_igblast_txt = rules.VDJB_igblast.output.VDJB_igblast_txt
        output: 
            VDJB_changeo_db = PVDJB + '/{sample}/changeo_db-pass.tsv',
            VDJB_changeo_clone = PVDJB + '/{sample}/changeo_clone-pass.tsv',
            VDJB_changeo_fail = PVDJB + '/{sample}/changeo_clone-fail.tsv',
            VDJB_changeo  = PVDJB + '/{sample}/changeo_clone-pass_germ-pass.tsv',
        log: e = Plog + '/VDJB_changeo/{sample}.e', o = Plog + '/VDJB_changeo/{sample}.o'
        benchmark: Plog + '/VDJB_changeo/{sample}.bmk'
        resources: cpus=config['VDJB_changeo_cpus']
        conda: 'envs/VDJ.yaml'
        params: 
            changeo_VB_ref = f'src/igblast-ref/{species}-VDJB/IGV.fasta',
            changeo_DB_ref = f'src/igblast-ref/{species}-VDJB/IGD.fasta',
            changeo_JB_ref = f'src/igblast-ref/{species}-VDJB/IGJ.fasta',
            outdir = PVDJB + '/{sample}'
        shell:"""
            MakeDb.py igblast -i {input.VDJB_igblast_txt} \\
                -r {params.changeo_VB_ref} {params.changeo_DB_ref} {params.changeo_JB_ref} \\
                -s {input.VDJB_orf_nt_fa} --outname changeo --regions default \\
                --failed --partial --format airr --extended 1>>{log.o} 2>>{log.e}
            
            DefineClones.py -d {output.VDJB_changeo_db} --failed --act set --nproc {resources.cpus}\\
                --outname changeo --model ham --norm len --dist 0.15 1>>{log.o} 2>>{log.e}

            # in case all contigs are failed or all passed
            touch {output.VDJB_changeo_clone}
            touch {output.VDJB_changeo_fail}

            CreateGermlines.py -d {output.VDJB_changeo_clone} -g dmask --cloned \\
                -r {params.changeo_VB_ref} {params.changeo_DB_ref} {params.changeo_JB_ref} \\
                1>>{log.o} 2>>{log.e}
            """

    rule VDJB_anarci:
        input: 
            VDJB_igblast_tsv = rules.VDJB_igblast.output.VDJB_igblast_tsv
        output:
            VDJB_orf_aa_fa = PVDJB + '/{sample}/seq_orf_aa.fasta',
            VDJB_gm_aa_fa = PVDJB + '/{sample}/seq_gm_aa.fasta',
            VDJB_orf_anarci_H = PVDJB + '/{sample}/orf_anarci_H.csv',
            VDJB_orf_anarci_KL = PVDJB + '/{sample}/orf_anarci_KL.csv',
            VDJB_gm_anarci_H = PVDJB + '/{sample}/gm_anarci_H.csv',
            VDJB_gm_anarci_KL = PVDJB + '/{sample}/gm_anarci_KL.csv',
            VDJB_ext_anarci_H = PVDJB + '/{sample}/ext_anarci_H.csv',
            VDJB_ext_anarci_KL = PVDJB + '/{sample}/ext_anarci_KL.csv',
        log: e = Plog + '/VDJB_anarci/{sample}.e', o = Plog + '/VDJB_anarci/{sample}.o'
        benchmark: Plog + '/VDJB_anarci/{sample}.bmk'
        resources: cpus=config['VDJB_anarci_cpus']
        params: 
            outdir = PVDJB + '/{sample}', 
            ext_numbering = config['ext_numbering'],
            VDJB_airr = rules.count.output.count_dir + f"/{config['count_VDJB_airr']}"
        conda: 'envs/VDJ.yaml'
        shell:"""
            # fetch sequences in ORF
            R -e " \
                x <- readr::read_tsv('{params.VDJB_airr}'); \
                genogamesh::parse_CellRanger_vdjseq(x, file='{output.VDJB_orf_aa_fa}', fa_content='seq_orf_aa'); \
                pre_y <- dplyr::pull(readr::read_tsv('{input.VDJB_igblast_tsv}'), germline_alignment_aa, sequence_id); \
                y <- stringr::str_replace_all(pre_y, '-', ''); names(y) <- names(pre_y); \
                genogamesh::write_fasta(y, '{output.VDJB_gm_aa_fa}'); \
                " 1>>{log.o} 2>>{log.e}

            cd {params.outdir}
            # for mutation target profile 
            ANARCI -i {output.VDJB_orf_aa_fa} -o orf_anarci -ht orf_anarci_hittable.txt \\
                --use_species {species} --restrict ig -s imgt --csv --ncpu {resources.cpus} \\
                --assign_germline 1>>{log.o} 2>>{log.e}

            ANARCI -i {output.VDJB_gm_aa_fa} -o gm_anarci -ht gm_anarci_hittable.txt \\
                --use_species {species} --restrict ig -s imgt --csv --ncpu {resources.cpus} \\
                --assign_germline 1>>{log.o} 2>>{log.e}

            # for extended numbering
            ANARCI -i {output.VDJB_orf_aa_fa} -o ext_anarci -ht ext_anarci_hittable.txt \\
                --use_species {species} --restrict ig -s {params.ext_numbering} --csv --ncpu {resources.cpus} \\
                --assign_germline 1>>{log.o} 2>>{log.e}
        """

    rule VDJB_parse:
        input:
            count_dir = rules.count.output.count_dir,
            VDJB_igblast_tsv = rules.VDJB_igblast.output.VDJB_igblast_tsv,
            VDJB_igblast_txt = rules.VDJB_igblast.output.VDJB_igblast_txt,
            VDJB_changeo = rules.VDJB_changeo.output.VDJB_changeo,
            VDJB_changeo_fail = rules.VDJB_changeo.output.VDJB_changeo_fail,
            VDJB_orf_anarci_H = rules.VDJB_anarci.output.VDJB_orf_anarci_H,
            VDJB_orf_anarci_KL = rules.VDJB_anarci.output.VDJB_orf_anarci_KL,
            VDJB_gm_anarci_H = rules.VDJB_anarci.output.VDJB_gm_anarci_H,
            VDJB_gm_anarci_KL = rules.VDJB_anarci.output.VDJB_gm_anarci_KL,
            VDJB_ext_anarci_H = rules.VDJB_anarci.output.VDJB_ext_anarci_H,
            VDJB_ext_anarci_KL = rules.VDJB_anarci.output.VDJB_ext_anarci_KL,
        output:
            VDJB_csv = PVDJB + '/{sample}/VDJB.csv',
            VDJB_mut = PVDJB + '/{sample}/VDJB_mut.csv',
            VDJB_stat = PVDJB + '/{sample}/VDJB_stat.yaml',
            stat_dir = directory(Pstat + '/{sample}/VDJB')
        log: notebook = Plog + '/VDJB_parse/{sample}.r.ipynb', 
             e = Plog + '/VDJB_parse/{sample}.e', o = Plog + '/VDJB_parse/{sample}.o'
        benchmark: Plog + '/VDJB_parse/{sample}.bmk'
        resources: cpus=config['VDJB_parse_cpus']
        conda: 'envs/VDJ.yaml'
        notebook: Plog + '/VDJB_parse.r.ipynb'


##################################
### VDJT
##################################
if lambda wildcards:config[wildcards.sample]['VDJT']:
    rule VDJT_parse:
        input:
            count_dir = rules.count.output.count_dir,
        output:
            VDJT_csv = PVDJT + '/{sample}/VDJT.csv',
            VDJT_stat = PVDJT + '/{sample}/VDJT_stat.yaml',
            stat_dir = directory(Pstat + '/{sample}/VDJT')
        log: notebook = Plog + '/VDJT_parse/{sample}.r.ipynb', 
             e = Plog + '/VDJT_parse/{sample}.e', o = Plog + '/VDJT_parse/{sample}.o'
        benchmark: Plog + '/VDJT_parse/{sample}.bmk'
        resources: cpus=config['VDJT_parse_cpus']
        conda: 'envs/VDJ.yaml'
        notebook: Plog + '/VDJT_parse.r.ipynb'



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
    input: filter_input = get_filter_input,
    output:
        filter_dir = directory(Pfilter + '/{sample}'),
        filter_stat = Pfilter + '/{sample}/filter_stat.yaml'
    params: 
        stat_dir = Pstat + '/{sample}',
        Bcell_changeo_flt_H = Pfilter + '/{sample}/Bcell_changeo_flt_H.tsv',
        Bcell_changeo_flt_L = Pfilter + '/{sample}/Bcell_changeo_flt_L.tsv'
    log: notebook = Plog + '/filter/{sample}.r.ipynb', 
         e = Plog + '/filter/{sample}.e', o = Plog + '/filter/{sample}.o'
    benchmark: Plog + '/filter/{sample}.bmk'
    resources: cpus=config['filter_cpus']
    conda: 'envs/visualize.yaml'
    notebook: Plog + '/filter.r.ipynb'


rule concat_stat:
    input: filter_input = get_filter_input
    output: stats = Pstat + '/{sample}/stats.csv'
    log: notebook = Plog + '/concat_stat/{sample}.r.ipynb', 
         e = Plog + '/concat_stat/{sample}.e', o = Plog + '/concat_stat/{sample}.o'
    params: metadata = config['metadata']
    benchmark: Plog + '/concat_stat/{sample}.bmk'
    resources: cpus=config['concat_stat_cpus']
    conda: 'envs/visualize.yaml'
    notebook: Plog + '/concat_stat.r.ipynb'

if config['run_batch_stat']:
    rule batch_stat:
        input: 
            Lstats = expand(Pstat + '/{sample}/stats.csv', sample=Lsample),
        output: 
            batch_stat = Pstat + '/batch_stat.csv',
        log: e = Plog + '/batch_stat.e', o = Plog + '/batch_stat.o'
        resources: cpus=1
        run: 
            pd.concat([pd.read_csv(str(f)) for f in input.Lstats], axis=0).to_csv(output.batch_stat, index=False)


##################################
### visualize
##################################
rule visualize:
    input: filter_dir = rules.filter.output.filter_dir
    output: visualize_rds = Pvisualize + '/{sample}.rds'
    log: notebook = Plog + '/visualize/{sample}.r.ipynb', 
         e = Plog + '/visualize/{sample}.e', o = Plog + '/visualize/{sample}.o'
    params: 
        stat_dir = Pstat + '/{sample}',
        echarts_theme = 'src/echarts_theme/mytheme.json',
        metadata = config['metadata']
    benchmark: Plog + '/visualize/{sample}.bmk'
    resources: cpus=config['visualize_cpus']
    conda: 'envs/visualize.yaml'
    notebook: Plog + '/visualize.r.ipynb'

rule visualize_rmd:
    input: visualize_rds = rules.visualize.output.visualize_rds
    output: visualize_html = Pvisualize + '/{sample}.html'
    log: e = Plog + '/visualize_rmd/{sample}.e', o = Plog + '/visualize_rmd/{sample}.o'
    params: css_dir = 'src/visualize_css'
    benchmark: Plog + '/visualize_rmd/{sample}.bmk'
    resources: cpus=config['visualize_rmd_cpus']
    conda: 'envs/visualize.yaml'
    script: 'scripts/visualize.Rmd'



##################################
### Bcell tree
##################################
if lambda wildcards: (config[wildcards.sample]['VDJB'] and config['run_tree']) :
    rule Bcell_tree:
        input:
            filter_dir = rules.filter.output.filter_dir
        output:
            Bcell_tree_H = Ptree + '/{sample}/Bcell_tree_H_igphyml-pass.tab',
            Bcell_tree_L = Ptree + '/{sample}/Bcell_tree_L_igphyml-pass.tab',
        log: e = Plog + '/Bcell_tree/{sample}.e', o = Plog + '/Bcell_tree/{sample}.o'
        benchmark: Plog + '/Bcell_tree/{sample}.bmk'
        resources: cpus=config['Bcell_tree_cpus']
        conda: 'envs/VDJ.yaml'
        params: 
            Bcell_changeo_flt_H = rules.filter.params.Bcell_changeo_flt_H,
            Bcell_changeo_flt_L = rules.filter.params.Bcell_changeo_flt_L,
            outdir = Ptree + '/{sample}',
            tree_sample = config['tree_sample']
        shell:"""
            BuildTrees.py -d {params.Bcell_changeo_flt_H} --collapse \\
                --sample {params.tree_sample} --igphyml --clean all --nproc {resources.cpus} \\
                --outdir {params.outdir} --outname Bcell_tree_H \\
                1>>{log.o} 2>>{log.e}
            BuildTrees.py -d {params.Bcell_changeo_flt_L} --collapse \\
                --sample {params.tree_sample} --igphyml --clean all --nproc {resources.cpus} \\
                --outdir {params.outdir} --outname Bcell_tree_L \\
                1>>{log.o} 2>>{log.e}
        """
