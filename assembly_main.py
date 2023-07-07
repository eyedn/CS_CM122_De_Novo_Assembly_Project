###############################################################################
#   Aydin Karatas
#   CS CM122 Project 3
#   project3_main.py
###############################################################################
import assembly_functions as asfx
from sys import argv
from collections import Counter
import datetime

def main():
    file = argv[1]
    input("Prove a k value. Press Enter to continue...")
    k = int(input("Enter the value for k: "))
    # read in reads from reads file
    print(f'{datetime.datetime.now()}: reading in reads')
    my_reads = asfx.read_reads(file)
    # create kmers from all reads
    print(f'{datetime.datetime.now()}: generating read {k}-mers')
    my_kmers_dict = asfx.kmerize_reads(my_reads, k)
    count_distribution = Counter(my_kmers_dict.values())
    with open("distribution.txt", "w") as f:
        sorted_counts = sorted(count_distribution.items(), key=lambda x: x[0])
        for count, frequency in sorted_counts:
            line = f"{count}: {frequency}\n"
            f.write(line)
    input("Data has been written to distribution.txt. Press Enter to continue...")
    min_occ = int(input("Enter the value for min_occ: "))
    my_kmers = []
    print(f'{datetime.datetime.now()}: removing erroneous {k}-mers')
    for kmer in my_kmers_dict:
        if my_kmers_dict[kmer] >= min_occ:
            my_kmers.append(kmer)
    # build adjecency matrix of db graph
    print(f'{datetime.datetime.now()}: creating de bruijn graph adjacency matrix')
    my_adj_matrix = asfx.bd_adjacency(my_kmers)
    # create a list of all outedges
    my_out_edges = asfx.all_out_edges(my_adj_matrix)
    # traverse the adjecency matrix to find a sequence
    print(f'{datetime.datetime.now()}: finding start nodes')
    inital_node = asfx.find_first_node(my_adj_matrix, my_out_edges)
    contigs = []
    print(f'{datetime.datetime.now()}: finding contigs with eulerian walks')
    for node in inital_node:
        my_eulerian_path, my_adj_matrix, contigs_IGNORE = asfx.traverse_adjacency(my_adj_matrix, node)
        contigs.append(asfx.trivial_string_spelling(my_eulerian_path))
    print(f'{datetime.datetime.now()}: pasting contigs together')
    de_nove_genome = ''.join(str(i) for i in contigs)
    with open("de_nove_genome.fasta", "w") as g:
        g.write(">de novo genome\n")
        g.write(de_nove_genome)

if __name__ == "__main__":
    main()