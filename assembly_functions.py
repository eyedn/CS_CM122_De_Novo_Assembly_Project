###############################################################################
#   Aydin Karatas
#   CS CM122 Project 3
#   project3_functions.py
###############################################################################
import assembly_classes as ascl
from typing import TextIO

# read in reads
def read_reads(file: TextIO) -> list[ascl.Read]:
    all_reads: list[ascl.Read] = []
    with open(file, "r") as f:
        for line in f:
            curr_line = line.strip()
            # ">" indicates the start of a new read
            if curr_line.startswith(">"):
                all_reads.append(ascl.Read(curr_line[1:]))  # Create a new ascl.Read instance
            else:
                # While we do not encounter a new ">", add the sequence to the most recent read
                all_reads[-1].add_to_sequence(curr_line)
    return all_reads

# trivially find the adjacency list of a de Bruij Graph
def bd_adjacency(kmers: list[str]) -> dict:
    adjacency = {}
    for kmer in kmers:
        edge = [kmer[1:], False]
        try:
            adjacency[kmer[:-1]].append(edge)
        except:
            adjacency[kmer[:-1]] = [edge]
    return adjacency

# def create kmer dict
def kmerize_reads(reads: list[ascl.Read], k: int) -> dict[str, int]:
    kmers:  dict[str, int] = {}
    for read in reads:
        read_kmers = create_kmers(read.sequence, k)
        for kmer in read_kmers:
            try:
                kmers[kmer] += 1
            except:
                kmers[kmer] = 1
    return kmers

# create kmers of a string
def create_kmers(text: str, k: int) -> list[str]:
    kmers = []
    for i in range(len(text)-k+1):
        kmers.append(text[i:i+k])
    return kmers

# traverse a full graph given the adjacency matrix
def traverse_adjacency(adjacency_matrix: dict, key_node = "") -> list[str]:
    contigs = []
    sequence = []
    nodes_with_out_edges = adjacency_matrix.keys()
    # initialize path with the first node
    i = 0 
    if key_node == None or key_node == "":
        key_node = nodes_with_out_edges[i]
    partial_sequence, adjacency_matrix, contigs = traverse_partial_adjacency(adjacency_matrix, key_node, 
                                                                    [key_node], contigs)
    sequence = partial_sequence
    while i < len(sequence):
        key_node = sequence[i]
        # if this node does not have any out_edges, it must be the laste node
        if key_node not in nodes_with_out_edges:
            # move to the next node 
            i += 1 
            continue
        adj_nodes = adjacency_matrix[key_node]
        # check for nodes in the sequence that what unvisted edges
        if any(adj_node[1] == False for adj_node in adj_nodes):
            # find the partial path of the unvisited edge
            partial_sequence, adjacency_matrix, contigs = traverse_partial_adjacency(adjacency_matrix, key_node, 
                                                                            [key_node], contigs)
            for j, s in enumerate(sequence):
                # location where to insert the new partial node is when a node
                # in the sequence equals the first node of the partial sequence
                if s != partial_sequence[0]:
                    continue
                sequence = sequence[:j] + partial_sequence + sequence[j+1:]
                break
            # stay on same index of sequence
        else:
            i += 1  # Move to the next key_node
    return(sequence, adjacency_matrix, contigs)

# traverse a subsection of the graph
def traverse_partial_adjacency(adjacency_matrix: dict, key: str, 
                               sequence: list[str], contigs: list[str]) -> str:
    contigs.append([key])
    first = True
    while True:
        if key not in adjacency_matrix.keys():
            contigs[-1].append(key)
            new_contig = trivial_string_spelling(contigs[-1])
            contigs[-1] = new_contig
            return sequence, adjacency_matrix, contigs
        outgoing_edges = adjacency_matrix[key]
        found_unvisited = False
        if len(outgoing_edges) > 1 and not first:
            contigs[-1].append(key)
            new_contig = trivial_string_spelling(contigs[-1])
            contigs[-1] = new_contig
            contigs.append([key])
            first = True
        for adj_node in outgoing_edges:
            if not adj_node[1]:
                adj_node[1] = True
                sequence.append(adj_node[0])
                if first:
                    first = False
                else:
                    contigs[-1].append(key)
                key = adj_node[0]
                found_unvisited = True
                break
        if not found_unvisited:
            new_contig = trivial_string_spelling(contigs[-1])
            contigs[-1] = new_contig
            return sequence, adjacency_matrix, contigs

# find starting node depending on in_edges - out_edges = -1
def find_first_node(adjecency_matrix: dict, out_edges: list[str]) -> list[str]:
    start_nodes = []
    for node in adjecency_matrix:
        out_degree = len(adjecency_matrix[node])
        in_degree = sum(edge == node for edge in out_edges)
        # return at the first instance of a semibalanced node
        if out_degree > in_degree:
            start_nodes.append(node)
    return start_nodes

# generate a list of all out-edges
def all_out_edges(adjecency_matrix: dict) -> list[str]:
    out_edges = []
    for node_edge in adjecency_matrix.values():
        for indiv_edge in node_edge:
            out_edges.append(indiv_edge[0])
    return(out_edges)

# trivially combine sequential kmers into 1 sting
def trivial_string_spelling(kmers: list):
    text = kmers[0]
    first = True
    for kmer in kmers:
        if first:
            first = False
            continue
        text += str(kmer)[-1]
    return text