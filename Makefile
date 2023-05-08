setup:
	bash init.sh && \
	mkdir -p ~/.config/snakemake/cell10x && \
	cp smk_profiles/cell10x.yaml ~/.config/snakemake/cell10x/config.yaml && \
	mkdir -p ~/.config/snakemake/cell10x_slurm && \
	cp smk_profiles/cell10x_slurm.yaml ~/.config/snakemake/cell10x_slurm/config.yaml && \
	snakemake --profile cell10x --conda-create-envs-only -j1
build_envs:
	snakemake --profile cell10x --conda-create-envs-only -j1
clean_envs:
	snakemake --profile cell10x --conda-cleanup-envs -j1
touch_all:
	snakemake --profile cell10x --touch -F -j1
touch_code_change:
	snakemake -R $(snakemake --list-code-changes) --touch