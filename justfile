build_envs:
	snakemake --profile profiles/local --conda-create-envs-only -j1
clean_envs:
	snakemake --profile profiles/local --conda-cleanup-envs -j1
touch_all:
	snakemake --profile profiles/local --touch -F -j1
touch_code_change:
	snakemake -R $(snakemake --list-code-changes) --touch
