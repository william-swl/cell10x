# in germline_db/human-VDJB

# igblast germline database
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKV.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLV.fasta
cat IGHV.fasta IGLV.fasta IGKV.fasta > IGV.fasta
rm IGHV.fasta IGLV.fasta IGKV.fasta

~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl IGV.fasta > human_V


# igblast blast database
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in human_V
cp human_V.* ../../internal_data/human/

# manual create igblast germline V gene annotation file
# human.ndm.imgt for igblastn, imgt domain_system
