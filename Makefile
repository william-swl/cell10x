setup:
	bash init.sh && \
	mkdir -p ~/.config/snakemake/cell10x && \
	cp smk_profiles/cell10x.yaml ~/.config/snakemake/cell10x/config.yaml && \
	snakemake --profile cell10x --conda-create-envs-only --until qc -j1
clean_envs:
	snakemake --profile cell10x --conda-cleanup-envs -j1
