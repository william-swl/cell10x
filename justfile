init:
	bash init.sh
	just build_envs
build_envs:
	conda config --set channel_priority strict
	snakemake --profile profiles/local --conda-create-envs-only -j1
clean_envs:
	snakemake --profile profiles/local --conda-cleanup-envs -j1
clean_conda_cache:
	conda clean --all
touch_all:
	snakemake --profile profiles/local --touch -F -j1
touch_code_change:
	snakemake -R $(snakemake --list-code-changes) --touch
