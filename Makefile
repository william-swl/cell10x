build_envs:
	snakemake --profile cell10x --conda-create-envs-only -j1
clean_envs:
	snakemake --profile cell10x --conda-cleanup-envs -j1
touch_all:
	snakemake --profile cell10x --touch -F -j1
touch_code_change:
	snakemake -R $(snakemake --list-code-changes) --touch