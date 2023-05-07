#!/usr/bin/bash
srcdir=src
log=$CONDA_PREFIX/deploy.log


##################################
### Install cellranger-7.1.0
##################################
cr_url="https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1683406556&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2ODM0MDY1NTZ9fX1dfQ__&Signature=iACtsRdEsNKtdOQwjCyLXXjn~1Kr0DmQii8XLHz8N~SRMpaSC3Sm7dSs5D3ajtyq-vBPbyw0Qa~TuYbmqMynFFao749Zujap~tNGwUquQzXjRXfgY-z7vPjlf2m~nD93zb5avPNjzXSuHWWYp1dbUt3r7ltqxiu9Khq4ZPnmLNgfKZMjDP6xrgqgXZCEIkpsPG1N49dPJqRjHas1IUR0IKfFJv2i6feInXbdpKpKmu2bcCNL6m9ygHvXybz4Yzx-RFDMIcjy3G71ADG86yI2Xa3RozEHDZSkRCNvIWi9i~HmcqPcQpRGQCbPVbC78SiJ15T635ifJn4QDhtO7NngHw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
cr_md5=0938af789631800b20176527063c029a

# installed
if [ -d "${srcdir}/cellranger-7.1.0" ]; then
    echo ">>> Cellranger already installed" >> $log 2>&1
else
    # downloaded
    if [ -f "${srcdir}/cellranger-7.1.0.tar.gz" ]; then
        md5_local=`md5sum ${srcdir}/cellranger-7.1.0.tar.gz | awk '{print $1}'`
        # download unfinished
        if [ "$md5_local" == "$cr_md5" ]; then
            echo ">>> Cellranger already downloaded" >> $log 2>&1
        # download finished
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


#############################################
### Download human gene expression reference
#############################################
gex_ref_url=https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
gex_ref_md5=886eeddde8731ffb58552d0bb81f533d
if [ -d "${srcdir}/human-gex-ref" ]; then
    echo ">>> Cellranger human gene expression reference already installed" >> $log 2>&1
else
    # downloaded
    if [ -f "${srcdir}/human-gex-ref.tar.gz" ]; then
        md5_local=`md5sum ${srcdir}/human-gex-ref.tar.gz | awk '{print $1}'`
        # download unfinished
        if [ "$md5_local" == "$gex_ref_md5" ]; then
            echo ">>> Cellranger human gene expression reference already downloaded" >> $log 2>&1
        # download finished
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
    echo ">>> Cellranger human vdj reference already installed" >> $log 2>&1
else
    # downloaded
    if [ -f "${srcdir}/human-vdj-ref.tar.gz" ]; then
        md5_local=`md5sum ${srcdir}/human-vdj-ref.tar.gz | awk '{print $1}'`
        # download unfinished
        if [ "$md5_local" == "$vdj_ref_md5" ]; then
            echo ">>> Cellranger human vdj reference already downloaded" >> $log 2>&1
        # download finished
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
