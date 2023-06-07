# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRAV.fasta
# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRAJ.fasta
# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBV.fasta
# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBD.fasta
# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBJ.fasta
# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRDV.fasta
# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRDD.fasta
# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRDJ.fasta
# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRGV.fasta
# wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRGJ.fasta

## delete duplicated TRAV/DV in TRDV.fasta

cat TRAV.fasta TRBV.fasta TRDV.fasta TRGV.fasta > TRV.fasta
rm TRAV.fasta TRBV.fasta TRDV.fasta TRGV.fasta
~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl TRV.fasta > human_gl_V
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in human_gl_V

cat TRBD.fasta TRDD.fasta > TRD.fasta
rm TRBD.fasta TRDD.fasta
~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl TRD.fasta > human_gl_D
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in human_gl_D


cat TRAJ.fasta TRBJ.fasta TRDJ.fasta TRGJ.fasta > TRJ.fasta
rm TRAJ.fasta TRBJ.fasta TRDJ.fasta TRGJ.fasta
~/software/ncbi-igblast-1.17.1/bin/edit_imgt_file.pl TRJ.fasta > human_gl_J
~/software/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in human_gl_J

