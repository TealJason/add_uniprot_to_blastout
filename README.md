This script 
1.runs blastp on a fasta files
2.filters the lowest E-value for each match for each subject from the blast results as a new file
3.uses the uniprot api to create a mapping file mapping each refseq id to a uniprot id and name
4.appends the uniprot id and name to the filtered blast results as new columns after the refseq id
