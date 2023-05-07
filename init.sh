#!/usr/bin/bash
current_dir=$(pwd)
sed -i "s|/home/william/pipeline/cell10x|$current_dir|g" sample_config/test.yaml
