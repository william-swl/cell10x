R.home()

baizer::pkglib(tidyverse, baizer)

sample <- snakemake@wildcards[['sample']]

# replace default configs with sample configs
config <- replace_item(snakemake@config, snakemake@config[[sample]])

# write feature_ref.csv
if (config[[sample]][['FB']]) {
    id2seq <- config[[sample]][['id2seq']]
    TBfb_ref <- map2_dfr(id2seq, names(id2seq), 
                         ~list2df(.x, colnames='sequence') %>% r2c('id') %>% mutate(type=.y)) %>%
        mutate(
            name = str_glue('{type}_{id}'),
            read = config[['FB_read']],
            pattern = config[['FB_pattern']],
            feature_type = config[['FB_type']]
              ) %>%
        select(id, name, read, pattern, sequence, feature_type)

    TBfb_ref %>% write_excel_csv(snakemake@params[['feature_ref']])
}

# write config.csv

out <- file(snakemake@output[['cellranger_config']], 'w')

# reference
writeLines(str_glue("[gene-expression]\nreference,{config[['gex_ref']]}\nchemistry,SC5P-R2"), out)
writeLines(str_glue("[feature]\nreference,{snakemake@params[['feature_ref']]}"), out)
writeLines(str_glue("[vdj]\nreference,{config[['vdj_ref']]}"), out)

# library
writeLines('[libraries]\nfastq_id,fastqs,lanes,feature_types,subsample_rate', out)
if (config[[sample]][['mRNA']]) {
    writeLines(str_glue("{sample}-mRNA,{config[['indir']]}/{sample},,Gene Expression,,"), out)
}
if (config[[sample]][['FB']]) {
    writeLines(str_glue("{sample}-FB,{config[['indir']]}/{sample},,{config[['FB_type']]},,"), out)
}
if (config[[sample]][['VDJB']]) {
    writeLines(str_glue("{sample}-VDJB,{config[['indir']]}/{sample},,VDJ-B,,"), out)
}
if (config[[sample]][['VDJT']]) {
    writeLines(str_glue("{sample}-VDJT,{config[['indir']]}/{sample},,VDJ-T,,"), out)
}

close(out) 