#!/usr/bin/bash

# edit profiles
current_dir=$(pwd)
sed -i "s|/home/william/pipeline/cell10x|$current_dir|g" smk_profiles/cell10x.yaml
sed -i "s|/home/william/pipeline/cell10x|$current_dir|g" smk_profiles/cell10x_slurm.yaml

# copy profiles
mkdir -p ~/.config/snakemake/cell10x
cp smk_profiles/cell10x.yaml ~/.config/snakemake/cell10x/config.yaml
mkdir -p ~/.config/snakemake/cell10x_slurm
cp smk_profiles/cell10x_slurm.yaml ~/.config/snakemake/cell10x_slurm/config.yaml

##################################
### input
##################################
read -p "Your Cellranger download link? (anykey to ignore)" cr_url
read -p "Download human reference? (y/n): " yn_human
read -p "Download mouse reference? (y/n): " yn_mouse

##################################
### Install cellranger-7.1.0
##################################
srcdir=src
mkdir -p $srcdir
log=./cellranger.log

cr_md5=0938af789631800b20176527063c029a

# installed
if [ -d "${srcdir}/cellranger-7.1.0" ]; then
    echo ">>> Cellranger has been installed" >> $log 2>&1
else
    # downloaded
    if [ -f "${srcdir}/cellranger-7.1.0.tar.gz" ]; then
        md5_local=`md5sum ${srcdir}/cellranger-7.1.0.tar.gz | awk '{print $1}'`
        # download finished
        if [ "$md5_local" == "$cr_md5" ]; then
            echo ">>> Cellranger has been downloaded" >> $log 2>&1
        # download unfinished
        else
            wget -c -O ${srcdir}/cellranger-7.1.0.tar.gz $cr_url >> $log 2>&1
            echo ">>> Cellranger downloaded" >> $log 2>&1
        fi
    # not downloaded
    else
        wget -c -O ${srcdir}/cellranger-7.1.0.tar.gz $cr_url >> $log 2>&1
        echo ">>> Cellranger downloaded" >> $log 2>&1
    fi

    # decompress
    tar -xzf ${srcdir}/cellranger-7.1.0.tar.gz -C ${srcdir}/ >> $log 2>&1
    echo ">>> Cellranger decompressed" >> $log 2>&1
    rm ${srcdir}/cellranger-7.1.0.tar.gz >> $log 2>&1
fi



if [ $yn_human = y ]; then
    #############################################
    ### Download human gene expression reference
    #############################################
    gex_ref_url=https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
    gex_ref_md5=dfd654de39bff23917471e7fcc7a00cd
    if [ -d "${srcdir}/human-gex-ref" ]; then
        echo ">>> Cellranger human gene expression reference has been installed" >> $log 2>&1
    else
        # downloaded
        if [ -f "${srcdir}/human-gex-ref.tar.gz" ]; then
            md5_local=`md5sum ${srcdir}/human-gex-ref.tar.gz | awk '{print $1}'`
            # download finished
            if [ "$md5_local" == "$gex_ref_md5" ]; then
                echo ">>> Cellranger human gene expression reference has been downloaded" >> $log 2>&1
            # download unfinished
            else
                wget -c -O ${srcdir}/human-gex-ref.tar.gz $gex_ref_url >> $log 2>&1
                echo ">>> Cellranger human gene expression reference downloaded" >> $log 2>&1
            fi
        # not downloaded
        else
            wget -c -O ${srcdir}/human-gex-ref.tar.gz $gex_ref_url >> $log 2>&1
            echo ">>> Cellranger human gene expression reference downloaded" >> $log 2>&1
        fi

        # decompress
        mkdir -p ${srcdir}/human-gex-ref
        tar -xzf ${srcdir}/human-gex-ref.tar.gz -C ${srcdir}/human-gex-ref --strip-components 1 >> $log 2>&1
        echo ">>> Cellranger human gene expression reference decompressed" >> $log 2>&1
        rm ${srcdir}/human-gex-ref.tar.gz >> $log 2>&1
    fi

    #######################################
    ### Download human VDJ reference
    #######################################
    vdj_ref_url=https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz
    vdj_ref_md5=65b5b033723b07fc1bb5375e5645761c
    if [ -d "${srcdir}/human-vdj-ref" ]; then
        echo ">>> Cellranger human vdj reference has been installed" >> $log 2>&1
    else
        # downloaded
        if [ -f "${srcdir}/human-vdj-ref.tar.gz" ]; then
            md5_local=`md5sum ${srcdir}/human-vdj-ref.tar.gz | awk '{print $1}'`
            # download finished
            if [ "$md5_local" == "$vdj_ref_md5" ]; then
                echo ">>> Cellranger human vdj reference has been downloaded" >> $log 2>&1
            # download unfinished
            else
                wget -c -O ${srcdir}/human-vdj-ref.tar.gz $vdj_ref_url >> $log 2>&1
                echo ">>> Cellranger human vdj reference downloaded" >> $log 2>&1
            fi
        # not downloaded
        else
            wget -c -O ${srcdir}/human-vdj-ref.tar.gz $vdj_ref_url >> $log 2>&1
            echo ">>> Cellranger human vdj reference downloaded" >> $log 2>&1
        fi

        # decompress
        mkdir -p ${srcdir}/human-vdj-ref
        tar -xzf ${srcdir}/human-vdj-ref.tar.gz -C ${srcdir}/human-vdj-ref --strip-components 1 >> $log 2>&1
        echo ">>> Cellranger human vdj reference decompressed" >> $log 2>&1
        rm ${srcdir}/human-vdj-ref.tar.gz >> $log 2>&1
    fi
fi



if [ $yn_mouse = y ]; then
    #############################################
    ### Download mouse gene expression reference
    #############################################
    gex_ref_url=https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
    gex_ref_md5=886eeddde8731ffb58552d0bb81f533d
    if [ -d "${srcdir}/mouse-gex-ref" ]; then
        echo ">>> Cellranger mouse gene expression reference has been installed" >> $log 2>&1
    else
        # downloaded
        if [ -f "${srcdir}/mouse-gex-ref.tar.gz" ]; then
            md5_local=`md5sum ${srcdir}/mouse-gex-ref.tar.gz | awk '{print $1}'`
            # download finished
            if [ "$md5_local" == "$gex_ref_md5" ]; then
                echo ">>> Cellranger mouse gene expression reference has been downloaded" >> $log 2>&1
            # download unfinished
            else
                wget -c -O ${srcdir}/mouse-gex-ref.tar.gz $gex_ref_url >> $log 2>&1
                echo ">>> Cellranger mouse gene expression reference downloaded" >> $log 2>&1
            fi
        # not downloaded
        else
            wget -c -O ${srcdir}/mouse-gex-ref.tar.gz $gex_ref_url >> $log 2>&1
            echo ">>> Cellranger mouse gene expression reference downloaded" >> $log 2>&1
        fi

        # decompress
        mkdir -p ${srcdir}/mouse-gex-ref
        tar -xzf ${srcdir}/mouse-gex-ref.tar.gz -C ${srcdir}/mouse-gex-ref --strip-components 1 >> $log 2>&1
        echo ">>> Cellranger mouse gene expression reference decompressed" >> $log 2>&1
        rm ${srcdir}/mouse-gex-ref.tar.gz >> $log 2>&1
    fi
    
    #######################################
    ### Download mouse VDJ reference
    #######################################
    vdj_ref_url=https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz
    vdj_ref_md5=c6f41db8f67aa83b04d64ba1ae96e681
    if [ -d "${srcdir}/mouse-vdj-ref" ]; then
        echo ">>> Cellranger mouse vdj reference has been installed" >> $log 2>&1
    else
        # downloaded
        if [ -f "${srcdir}/mouse-vdj-ref.tar.gz" ]; then
            md5_local=`md5sum ${srcdir}/mouse-vdj-ref.tar.gz | awk '{print $1}'`
            # download finished
            if [ "$md5_local" == "$vdj_ref_md5" ]; then
                echo ">>> Cellranger mouse vdj reference has been downloaded" >> $log 2>&1
            # download unfinished
            else
                wget -c -O ${srcdir}/mouse-vdj-ref.tar.gz $vdj_ref_url >> $log 2>&1
                echo ">>> Cellranger mouse vdj reference downloaded" >> $log 2>&1
            fi
        # not downloaded
        else
            wget -c -O ${srcdir}/mouse-vdj-ref.tar.gz $vdj_ref_url >> $log 2>&1
            echo ">>> Cellranger mouse vdj reference downloaded" >> $log 2>&1
        fi

        # decompress
        mkdir -p ${srcdir}/mouse-vdj-ref
        tar -xzf ${srcdir}/mouse-vdj-ref.tar.gz -C ${srcdir}/mouse-vdj-ref --strip-components 1 >> $log 2>&1
        echo ">>> Cellranger mouse vdj reference decompressed" >> $log 2>&1
        rm ${srcdir}/mouse-vdj-ref.tar.gz >> $log 2>&1
    fi

fi



##################################
### build envs
##################################
snakemake --profile cell10x --conda-create-envs-only -j1