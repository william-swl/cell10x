wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGHV.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGKV.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGLV.fasta
cat IGHV.fasta IGLV.fasta IGKV.fasta > IGV.fasta
rm IGHV.fasta IGLV.fasta IGKV.fasta
~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl IGV.fasta > mouse_gl_V
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in mouse_gl_V

wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGHD.fasta
mv IGHD.fasta IGD.fasta
~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl IGD.fasta > mouse_gl_D
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in mouse_gl_D

wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGHJ.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGKJ.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGLJ.fasta
cat IGHJ.fasta IGLJ.fasta IGKJ.fasta > IGJ.fasta
rm IGHJ.fasta IGLJ.fasta IGKJ.fasta
~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl IGJ.fasta > mouse_gl_J
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in mouse_gl_J
