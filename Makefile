setup:
	snakemake --profile cell10x --conda-create-envs-only --until qc -j1
clean_envs:
	snakemake --profile cell10x --conda-cleanup-envs -j1