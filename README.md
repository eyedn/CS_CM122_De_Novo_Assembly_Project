# CS CM122 Assembly Project
I created an assembly program that performs de novo genome assembly of reads data using an adjenct list to represent a Bruijn Graph. The project utilizes code from the Read Mapper Project to produce the output. 

The code for project 3a and 3b are the same. Here are the instructions.
1. Use the reads fasta file to create a de novo genome assembly. Use the command

python3 assembly_main.py <reads_file.fasta>

The code will promp the user to input a k value. I chose 20 for the highest score on the leaderboard I achieved. The code will then produce a file called <distributions.txt>, which contains the distribution for kmer occurances. The script will then prompt the user to input a minimum occurance threshold. I chose 12 for the highest score on the leaderboard. The file <de_novo_genome.fasta> will be created.

2.  Align the reads back to the <de_novo_genome.fasta> to get the order of reads. The code for this section was adapted from my burrows-wheeler implementation of project 1. Use the command

python3 alignment_main.py <reads_file.fasta>

This will produced the ordered list of reads in the file <predictions.csv>.
