wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKV.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLV.fasta
cat IGHV.fasta IGLV.fasta IGKV.fasta > IGV.fasta
rm IGHV.fasta IGLV.fasta IGKV.fasta
~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl IGV.fasta > human_gl_V
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in human_gl_V 

wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta
mv IGHD.fasta IGD.fasta
~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl IGD.fasta > human_gl_D
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in human_gl_D

wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHJ.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKJ.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLJ.fasta
cat IGHJ.fasta IGLJ.fasta IGKJ.fasta > IGJ.fasta
rm IGHJ.fasta IGLJ.fasta IGKJ.fasta
~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl IGJ.fasta > human_gl_J
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in human_gl_J
