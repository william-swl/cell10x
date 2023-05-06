# /usr/bin/bash
log=$CONDA_PREFIX/deploy.log

# Download cellranger
cr_url="https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1683406556&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2ODM0MDY1NTZ9fX1dfQ__&Signature=iACtsRdEsNKtdOQwjCyLXXjn~1Kr0DmQii8XLHz8N~SRMpaSC3Sm7dSs5D3ajtyq-vBPbyw0Qa~TuYbmqMynFFao749Zujap~tNGwUquQzXjRXfgY-z7vPjlf2m~nD93zb5avPNjzXSuHWWYp1dbUt3r7ltqxiu9Khq4ZPnmLNgfKZMjDP6xrgqgXZCEIkpsPG1N49dPJqRjHas1IUR0IKfFJv2i6feInXbdpKpKmu2bcCNL6m9ygHvXybz4Yzx-RFDMIcjy3G71ADG86yI2Xa3RozEHDZSkRCNvIWi9i~HmcqPcQpRGQCbPVbC78SiJ15T635ifJn4QDhtO7NngHw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
cr_md5=0938af789631800b20176527063c029a

# installed
if [ -d "src/cellranger-7.1.0" ]; then
    echo ">>> Cellranger already installed" >> $log 2>&1
else
    # downloaded
    if [ -f "src/cellranger-7.1.0.tar.gz" ]; then
        md5_local=`md5sum src/cellranger-7.1.0.tar.gz | awk '{print $1}'`
        # download unfinished
        if [ "$md5_local" == "$cr_md5" ]; then
            echo ">>> Cellranger already downloaded" >> $log 2>&1
        # download finished
        else
            wget -c -O src/cellranger-7.1.0.tar.gz $cr_url >> $log 2>&1
            echo ">>> Cellranger downloaded" >> $log 2>&1
        fi
    # not downloaded
    else
        wget -c -O src/cellranger-7.1.0.tar.gz $cr_url >> $log 2>&1
        echo ">>> Cellranger downloaded" >> $log 2>&1
    fi

    # decompress
    tar -xzf src/cellranger-7.1.0.tar.gz -C src/ >> $log 2>&1
    echo ">>> Cellranger decompressed" >> $log 2>&1
    rm src/cellranger-7.1.0.tar.gz >> $log 2>&1
fi
