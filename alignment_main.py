###############################################################################
#   Aydin Karatas
#   CS 122 S23
#   Project 1a
#   project1a_main.py
###############################################################################
import alignment_classes as alcl
import alignment_functions as alfx
from sys import argv
import datetime

def main():
    genome_file = "de_nove_genome.fasta"
    reads_file = argv[1]
    # create genome objcet and dictionary of reads objects
    print(f'{datetime.datetime.now()}: reading in genome')
    my_genome = alcl.Genome(genome_file)
    # reading in reads and mapping them to the genome
    my_reads, read_mutations_and_errors = alfx.create_reads_list(reads_file, 10, 7, my_genome, 3)
    # order reads by position
    my_ordered_reads = alfx.order_reads_by_mapped_idx(my_reads)
    with open("predictions.csv", "w") as f:
        for read in my_ordered_reads:
            f.write(read + '\n')


if __name__ == "__main__":
    main()