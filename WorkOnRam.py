 
#!/usr/bin/env python
# -*- coding: utf-8 -*-


def parse_args():
    '''Parse arguments given to script'''

    import argparse

    parser = argparse.ArgumentParser(description="Given two genomes, get both "\
        "common and specific regions to both genomes following a kmer-based algorithm")
    
    parser.add_argument("-d", dest="genome_dir", metavar="Genome dir", required=True, 
        help="Directory containing fasta or multifasta files to be clustered and a cinfo file")
    
    parser.add_argument("-r", dest="recursive", action="store_true", 
        help="Use for check all directories within directory provided by << -d >>")
    
    parser.add_argument("-g", metavar="jumping Gap between k-mers. DEFAULT: 6", dest="G", 
        default=6, help="Gap between k-mers. From inf to -k+1. (e.g G = 0 if no"\
        " gaps are desired. G = -K+1 if an overlap of 1 base between k-mers is desired")
    
    parser.add_argument("-k", metavar="k-mer length. DEFAULT: 12", dest="k", default=12, 
        help="k-mer length. Avoid use k > 12 in desktop computers")
    
    parser.add_argument("-f", metavar="min alignment length. DEFAULT: 120", dest="F", default=120,
        help="Distance checked by k-mer jumping from a seed k-mer to be considered as" \
        " a significant alignment")
    
    parser.add_argument("-j", metavar="max distance to combine fragments. DEFAULT: 25", 
        dest="J", default=25, help="Maximum distance to join significant alignments "\
        "from k-mer jumping")
    
    parser.add_argument("--max-seeds", metavar="maximum seeds per kmer. DEFAULT: 20", 
                        dest="max_seeds", default=20, help="Maximum seeds of a k-mer"\
                        " in order to seed alignments. This limit allows for a dramatic "\
                        "speedup in repetitive regions")
    
    parser.add_argument("-l, --length", dest="L", metavar="min cluster length. DEFAULT: 500", 
        default=500, help="Minimum length of a sequence in order to be considered as a new Cluster")
    
    parser.add_argument("-n", dest="N", metavar="max percentage of N per cluster. DEFAULT: 0.8",
     default=0.8, help="Maximum percentage of ambiguous nucleotides in clusters. If a new"\
     " cluster exceeds this maximum, it is discarded")

    parser.add_argument("-o", dest="out_orphan", action="store_true",
        help="Output orphan clusters that did not pass L and N filters"\
        "(.orphan and .omapping)")

    args = parser.parse_args()

    return args

def Align(k_start, seed_coordinate, index, sequence, k, G, reverse=False):
    '''After a seeding match, get kmers from sequence moving k+G bases between
    each one, and "moving" k+G bases from the seed_coordinate too. 
    For each n-k-mer, if it has a coordinate in index that coincides with 
    seed_coordinate + (k+G)*n, they are considered to be contiguous in both
    sequences, the scanned one and the indexed one. This function looks for all
    contiguous k-mers and returns the length of the extension produced

    Since each sequence has reverse coordinates for its reverse complement in
    index. Align also checks if contiguity is produced from the k_start and
    seed_coordinate backwards

    Example:
    A sequence with 12 bases with a k = 4 and G = 2 will have in index:

    For direct strand:
    AAAATTGGGGCC --> index = { "AAAA" : [0, ], "GGGG" : [6, ]

    And for reverse complement:
    GGCCCCAATTTT --> index = { "CCCC" : [6, ], [TTTT, 0]}
    '''
    from pang.seq_utils import GappedKmerGenerator
    import sys

    gapped_kmer_gen = GappedKmerGenerator(sequence, k_start, k, G)
    current_index_coord = seed_coordinate
    kmer = gapped_kmer_gen.next()
    jump = k + G
    # If alignment is being produced in reverse_complement, move "jump" backwards
    if reverse:
        jump = -jump
    while kmer in index:
        kmer_index_coords = index[kmer] # Take coords for next k+G gapped kmer
        if kmer_index_coords:
            current_index_coord += jump # Move k+G bases of index
            # When both coordinates don't coincide they are not contiguous
            if current_index_coord in kmer_index_coords:
                #print "current index_coord {} WAS in kmer_index_coords {}={}".format(current_index_coord, kmer, kmer_index_coords)

                # If both coordinates coincide, get and check next gapped kmer
                kmer = gapped_kmer_gen.next()
            else:
                #print "current index_coord {} WAS NOT in kmer_index_coords {}={}".format(current_index_coord, kmer, kmer_index_coords)
                current_index_coord -= jump
                if not reverse:
                    alignment_length = current_index_coord - seed_coordinate + k 
                else:
                    alignment_length = current_index_coord - seed_coordinate - k
                return alignment_length
        else:
            current_index_coord -= jump
            if not reverse:
                alignment_length = current_index_coord - seed_coordinate + k 
            else:
                alignment_length = current_index_coord - seed_coordinate - k
            
            return alignment_length
    
    # If kmer is returned as an int, an offset value is returned because the
    # end of sequence has been reached while scanning for kmers, check then
    # the last kmer
    if type(kmer) == int:
        last_offset = kmer
        last_right_coord = current_index_coord
        last_kmer = sequence[len(sequence)-k:]
        if last_kmer in index:
            kmer_index_coords = index[last_kmer]
            if kmer_index_coords:
                # Move the current index coord taking into account the offset
                current_index_coord += jump
                if not reverse:
                    current_index_coord -= last_offset
                else:
                    current_index_coord += last_offset
                if current_index_coord not in kmer_index_coords:
                    current_index_coord = last_right_coord
                if not reverse:
                    alignment_length = current_index_coord - seed_coordinate + k 
                else:
                    alignment_length = current_index_coord - seed_coordinate - k
                
                return alignment_length
    
    # If kmer not in index and its type is not an int, it is because
    # that kmer has ambiguous nucleotides, in that case return the alignment
    # made up to this point
    # 
    # This line is only reached if either, the kmer has ambiguous nucleotides or
    # it has no coords so kmer_index_coords = None
    if not reverse:
        alignment_length = current_index_coord - seed_coordinate + k 
    else:
        alignment_length = current_index_coord - seed_coordinate - k

    return alignment_length

def ExtendSeeds(k_start, seed_coordinates, index, sequence, k, G, F):
    '''This function takes all coordinates where a kmer of the scanned sequence
    is found in the indexed sequence and tries to make alignments starting from
    each seeding_coordinate. It returns a tuple of coordinates for the shortest
    aligned region of the scanned sequence and a list containing tuples of 
    coordinates for each aligment produced in the indexed sequence

    Parameter F allows to set a cutoff of a minimum extension for an aligment
    to be considered

    We keep the shortest alignment in order to continue scanning from the
    shortest fragment aligned to avoid pass some other k-mers that could seed
    another interesting alignments
    '''

    indexed_alignments = [] # To keep all aligments produced in indexed sequence
    shortest_alignment = int # int always > any value
    alignment = False # To check if any aligment > F has been produced
    alignments = [(-1,-1)] # To store coords for those regions already aligned

    # Two equal code blocks, first to align in forward index, second to align in
    # reverse_complement index

    # Align in forward index
    for seed_coordinate in seed_coordinates:
        if CheckSeed(seed_coordinate, alignments):
            alignment_length = Align(k_start, seed_coordinate, index, sequence, k, G)
            if alignment_length >= F:
                if alignment_length < shortest_alignment:
                    shortest_alignment = alignment_length # Keep shortest alignment
                indexed_start = seed_coordinate
                indexed_end = seed_coordinate + alignment_length - 1
                indexed_alignments.append( (indexed_start, indexed_end) )
                alignments.append((indexed_start, indexed_end))
                alignment = True

    # Align in reverse_complement index. Recall that alignment_length returned by
    # Align will be negative. e.g -160, -1250
    # ??? First reverse seed_coordinates list to avoid extension like okazaki fragments
    # ??? seed_coordinates.reverse()
    for seed_coordinate in seed_coordinates:
        if CheckSeed(seed_coordinate, alignments):
            alignment_length = Align(k_start, seed_coordinate, index, sequence, k, G,
                                     reverse=True)
            if -alignment_length >= F: # Recall that -(-160) is 160
                if -alignment_length < shortest_alignment:
                    shortest_alignment = -alignment_length # Keep shortest alignment
                indexed_start = seed_coordinate
                indexed_end = seed_coordinate + alignment_length + 1
                indexed_alignments.append( (indexed_start, indexed_end) )
                alignments.append((indexed_start, indexed_end))
                alignment = True

    if alignment:
        scanned_start = k_start
        scanned_end = k_start + shortest_alignment - 1
        scanned_alignment = (scanned_start, scanned_end)
        return (scanned_alignment, indexed_alignments)  
    else:
        return 0

def CheckSeed(seed, alignments):
    '''This function takes a seed coord and a list of tuples of coordinates
    with those alignments that have been already produced. It eliminates
    seeds that fall in regions already aligned:

    |----(Seed 1)----------(Seed 2)----------------------|
         |---------------------------------------------->| Alignment from Seed 1
         
         Then eliminate Seed 2 TO AVOID THIS:
                           |---------------------------->| Alignment from Seed 2
                                                           (shortest_alignment)
    
    This function returns False if Seed falls in a region already aligned
    or True otherwise
    '''

    for alignment_coords in alignments:
        start = alignment_coords[0]
        end = alignment_coords[1]
        if seed >= start and seed <= end:
            return False

    return True

def SeedAndExtend(sequence, index, k, G, F, max_seeds, k_start=0):

    '''This function takes k-mers from a sequence and looks if they can seed 
    alignments with the "indexed sequence". If so, it tries to extend those
    alignments. If alignments of length > F are produced, then coordinates of
    that alignments are added to the "core coordinates" and continues scanning
    from the end of the shortest_alignment produced in the scanned sequence

    It returns two lists for scanned and indexed sequence that contains tuples
    of pair of coordinates (start, end) that correspond to core regions in both
    sequences'''
    
    from pang.seq_utils import KmerGenerator, SkipAmbiguous

    alignment_coordinates = []
    non_ambiguous={"A", "T", "G", "C"}
    min_non_ambiguous = k # Parameter for SkipAmbiguous
    kmer_gen = KmerGenerator(sequence, k_start, k) # Get kmers overlaping by 1

    kmer = kmer_gen.next() # Get first kmer
    while True:
        if kmer in index:
            seed_coordinates = index[kmer] # Get coordinates of that kmer in index
            if seed_coordinates and len(seed_coordinates) <= max_seeds: 
                # Try to obtain aligments for each seed
                alignments = ExtendSeeds(k_start, seed_coordinates, index, sequence,
                                         k, G, F)
                if alignments:
                    alignment_coordinates.append(alignments)
                    scanned_alignment = alignments[0]
                    # k_start now is last coordinate of the shortest alignment
                    k_start = scanned_alignment[1] + 1
                    # Start again the generator with updated k_start
                    kmer_gen = KmerGenerator(sequence, k_start, k)
                    kmer = kmer_gen.next() # And get new kmer
                else: # If not alignment produced
                    kmer = kmer_gen.next()
                    k_start += 1
            else:
                kmer = kmer_gen.next() # If not seed coordinates
                k_start += 1
        else:
            # If kmer not in index, check if it is != 0, in that case
            # kmer is a string not in index, therefore is a kmer containing
            # ambiguous DNA bases. In that case, pass over the ambiguous sequence
            # until at least a number of contiguous nucleotides defined by 
            # min_non_ambiguous are found and continue checking from that point
            if kmer != 0:
                # Check position after contiguous non_ambiguous nucleotides
                k_start = SkipAmbiguous(sequence, k_start, non_ambiguous, 
                                        min_non_ambiguous)
                # Start again the generator with updated k_start
                kmer_gen = KmerGenerator(sequence, k_start, k)
                kmer = kmer_gen.next() # And get new kmer
            # If kmer is == 0 then the generator has finished and the end of
            # the sequence has been reached. In that case, return
            else:
                return alignment_coordinates

def NewMappingGroup(mapping_file, ref):
    '''This function updates the mapping_file after new sequences have ben added
    to the core genome file'''
    
    with open(mapping_file, "a") as handle:
        handle.write(str(CURRENT_CLUSTER) + "\t" + ref + "\n")

def UpdateMappingGroups(mapping_file, groups, ref):
    '''This function updates the mapping file by adding current ref to each group
    it has produced an alignment with'''
    import shutil

    mapping_file_tmp = mapping_file + ".tmp"
    with open(mapping_file) as handle:
        # Open a temporary file for update the current mapping_file info
        with open(mapping_file_tmp, "w") as tmp:
            for line in handle:
                group, refs = line.split("\t")
                # Write an exact line if that group is not implied
                if group not in groups:
                    tmp.write(line)
                else:
                    # If not, add ;ref to the end of the line
                    tmp.write(line.rstrip() + ";" + ref + "\n")
    # And mv the overwrite the older version with the updated one
    shutil.move(mapping_file_tmp, mapping_file)

def AlignRecords(fasta, index, k, G, F, J, L, N, max_seeds, pangenome, mapping, 
                  orphan_seqs, orphan_map, index_map):
    '''This function iterates over each record in the fasta file and performs
    the SeedAndExtend function over each one, returning coordinates of "core"
    alignments for each record and indexed sequence in a dict where record
    title is the key'''
    import sys
    from pang.parse_utils import FastaParser, GetID, GetRecordGroups
    from pang.coordinates_utils import SortCoordinates, JoinFragments
    from pang.coordinates_utils import  MapAlignments, GetNewCoreSeqs
    from pang.index_utils import ReindexRecord
    from pang.seq_utils import NucleotideFreq

    # NOTE that each time a coordinate is going to be appended to a mapping dict
    # we sum +1 so we switch from programming to human coordinates
    title = "FILE_EMPTY"
    cluster_title = pangenome["TITLE"]
    with open(fasta) as handle:
        for title, seq in FastaParser(handle):
            # Get the ref from current record
            ref = GetID(title)
            if seq:
                # First produce alignments with indexed sequences and retrieve
                # coordinates of these alignments
                alignment_coordinates = SeedAndExtend(seq, index, k, G, F, max_seeds)
                if alignment_coordinates:
                    # Sort alignment coordinates
                    sorted_coordinates = SortCoordinates(alignment_coordinates)
                    # And Join fragments that are as close as parameter J
                    joined_coords = JoinFragments(sorted_coordinates, J)
                    # Get new sequences that do not produce core alignments
                    # This new sequences will be added to a new GROUP so first
                    # update the global variable with the new group
                    new_seqs = GetNewCoreSeqs(joined_coords, seq, L, N)
                    for new_seq, new_seq_start, new_seq_end, orphan in new_seqs:
                        if not orphan:
                            curr_cluster = pangenome["CURRENT"]
                            cluster = "Cluster_" + str(curr_cluster)
                            # Format header to CORE_TITLE format
                            header = cluster_title + cluster
                            # Add new sequences to index and update index_map
                            index = ReindexRecord(header, k, index, index_map, new_seq)
                            # Update pangenome dictionary with new_core_seq
                            pangenome[header] = new_seq
                            # Update mapping_dict with new groups
                            ref_coords = (1, len(new_seq), ref, "plus", 
                                         new_seq_start + 1, new_seq_end + 1)
                            mapping[cluster] = [ref_coords]
                            pangenome["CURRENT"] += 1
                        else:
                            curr_orphan = orphan_seqs["CURRENT"]
                            curr_orphan = "Orphan_" + str(curr_orphan)
                            header = cluster_title + curr_orphan
                            orphan_seqs[header] = new_seq
                            orphan_coords = (1, len(new_seq), ref, "plus", 
                                             new_seq_start + 1, new_seq_end + 1)
                            orphan_map[curr_orphan] = [orphan_coords]
                            orphan_seqs["CURRENT"] += 1



                else:
                    # If no alignment is produced, scanned_core_coords is evaluated
                    # as False since it contains an empty list of coordinates
                    # in that case all sequence is new
                    if len(seq) > L and NucleotideFreq(seq, "N") < N:
                        curr_cluster = pangenome["CURRENT"]
                        cluster = "Cluster_" + str(curr_cluster)
                        header = cluster_title + cluster
                        index = ReindexRecord(header, k, index, index_map, seq)
                        # Update pangenome with a new_core_seq
                        pangenome[header] = seq
                        # Update mapping dict for that new_seq
                        ref_coords = (1, len(seq), ref, "plus", 1, len(seq))
                        mapping[cluster] = [ref_coords]
                        pangenome["CURRENT"] += 1
                        # If all sequence is new, there is no need to look which records
                        # each alignment maps to, so continue with following iteration of
                        # for loop
                        continue
                    else:
                        curr_orphan = orphan_seqs["CURRENT"]
                        curr_orphan = "Orphan_" + str(curr_orphan)
                        header = cluster_title + curr_orphan
                        orphan_seqs[header] = seq
                        orphan_coords = (1, len(seq), ref, "plus", 1, len(seq))
                        orphan_map[curr_orphan] = [orphan_coords]
                        orphan_seqs["CURRENT"] += 1
                        continue

            else: # If there is one empty seq
                sys.exit("Error: one or more records are empty")

            # For sequences that were already aligned in the core genome
            # check which clusters they have been aligned with
            mapped_alignments = MapAlignments(joined_coords, index_map)
            # If any core alignment has been produced
            if mapped_alignments:
                for scan_coords, clusters_mapped in mapped_alignments:
                    scan_start = scan_coords[0] + 1
                    scan_end = scan_coords[1] + 1
                    for cluster_mapped in clusters_mapped:
                        cluster = cluster_mapped[0]
                        c_start = cluster_mapped[1] + 1
                        c_end = cluster_mapped[2] + 1
                        if c_start < c_end:
                            strand = "plus"
                        else:
                            strand = "minus"
                        mapping[cluster].append(
                            (c_start, c_end, ref, strand, scan_start, scan_end )                            )

 
    if title == "FILE_EMPTY": # If no title, seq returned by parser title
                              # remains <<FILE_EMPTY>>
        sys.exit("Error: scanned fasta file seems to be empty")

    return pangenome, mapping, orphan_seqs, orphan_map

def InitCluster(cluster_info, genome_dir_path):
    '''This function initializes all files necessary to start build the core
    genome within a given directory. It gets the bigger fasta file as de first 
    indexed-reference file and creates a ".core" file containing that genome
    as the first record of the core genome belonging to Group 1.
    It also creates the .mapping file which will be initialize with Group 1 
    in first field and gi from that first reference in second field.

    ir returns .core and .mapping filenames with absolute path

    '''
    import os
    from pang.parse_utils import GetID, FastaParser
    
    genus = cluster_info[0]
    species = cluster_info[1]
    cluster_title = cluster_info[2]

    # Create the genus_species.core string to create/stat the file
    clusters_file = genus + "_" + species + ".clusters"
    mapping_file = genus + "_" + species + ".cmapping"
    orph_seq_file = genus + "_" + species + ".orphan"
    orph_map_file = genus + "_" + species + ".omapping"
    # Add absolute path to each file
    clusters_file = os.path.normpath(os.path.join(genome_dir_path, clusters_file))
    mapping_file = os.path.normpath(os.path.join(genome_dir_path, mapping_file))
    orph_seq_file = os.path.normpath(os.path.join(genome_dir_path, orph_seq_file))
    orph_map_file = os.path.normpath(os.path.join(genome_dir_path, orph_map_file))

    return clusters_file, mapping_file, orph_seq_file, orph_map_file, cluster_title

def Clusterize(genome_dir_path, k, G, F, J, L, N, index, max_seeds, out_orphan):
    '''This function takes a genome dir -that sould be a directory containing
    genome files in fasta format for a given specie- to build the set of 
    sequences that will form the core genome

    It creates:
    -a .mapping file that as a tab delimited text file containing a group number
    in the first field and all gi's pretaining to that group in the second field

    -a .core file, which is actually a fasta formated file containing all
    sequences that form the core genome

    Then it calls all necessary functions in order to align genomes and build
    the actual core

    It uses a core.info file that should have been created previously by the script
    that splits genomes from refseq .fna files in respective species folders.

    core.info has a row with common header for every sequence that will form the
    core, and another row with the number of different core groups
    '''

    from os.path import normpath, realpath
    from os.path import join as pjoin
    from pang.pang_IO import LoadInfo, WritePangenome, WriteMapping, ListFasta
    import sys
    import glob

    # Get full and normalized file paths
    genome_dir_path = realpath(genome_dir_path)
    # Load genus, species and title info from khan.info file
    info = normpath(pjoin(genome_dir_path, "cluster.info"))
    # If cluster.info exists
    if glob.glob(info):
        # Load its information
        cluster_info = LoadInfo(info)
        # Define output file names and cluster_title for the header of cluster seqs
        clusters_file, mapping_file, orph_seq_file, orph_map_file, cluster_title \
        = InitCluster(cluster_info, genome_dir_path)
        # Get al fasta files within the directory to be processed
        fasta_list = ListFasta(genome_dir_path)
        
        # Initialize dictionaries that will contain the actual clusters with
        # its sequences and auxiliar dictionaries with cluster alignment entries
        # to build mapping files in BED format
        # pangenome and orphan seqs are dictionaries that contain clusters ids as 
        # keys and sequences as values.
        # mapping and orphan_map contain entries for cluster alignments in tuples
        # as values and cluster ids as keys, so mapping files in BED format can be
        # written
        pangenome = {"CURRENT" : 1, "TITLE" : cluster_title}
        mapping = {}
        orphan_seqs = {"CURRENT" : 1, "TITLE" : cluster_title}
        orphan_map = {}
        # Index map contains both correlated lists of cluster ids and its respective
        # coordinates in index
        index_map = [ [] , [] ]
        # Perform alignments of records for each of the remaining genomes in list
        for fasta in fasta_list:
            print "Calculating pangenome for {}...".format(fasta)
            # Get absolute path to each genome
            fasta = normpath(pjoin(genome_dir_path, fasta))
            # And update pangenome and mapping dicts
            pangenome, mapping, orphan_seqs, orphan_map = \
            AlignRecords(fasta, index, k, G, F, J, L, N, max_seeds, pangenome, mapping,
                         orphan_seqs, orphan_map, index_map)
            # Once all records have been aligned, write final core and mapping from
            # pangenome and mapping dicts
        WritePangenome(pangenome, clusters_file)
        WriteMapping(mapping, mapping_file )

        if out_orphan:
            WritePangenome(orphan_seqs, orph_seq_file)
            WriteMapping(orphan_map, orph_map_file)
    else:
        sys.exit("Unable to find cluster.info file\n")



def ProcessDir(genome_dir, k, G, F, J, L, N, max_seeds, out_orphan):
    '''This function takes clustering parameters and a directory containing fasta files
    and calls the main function Cluster printing some info about processing time'''
    import os
    from time import time
    from pang.index_utils import BuildIndex

    genome_dir = os.path.realpath(genome_dir)
    start_time = time()
    index = BuildIndex(k)
    Clusterize(genome_dir, k, G, F, J, L, N, index, max_seeds, out_orphan)
    time_end = time() - start_time
    print "Pangenome calculated in {} seconds".format(time_end)

def ProcessDirRecursively(genomes_dir, k, G, F, J, L, N, max_seeds):
    '''This function takes a directory and performs the clustering in each of 
    the directories that it contains. It is meant to be used with a directory 
    that contains one folder per species when clustering whole databases so
    the index of size k has to be built only once and is reused for each species'''
    import os
    from time import time
    from pang.index_utils import BuildIndex
    
    genomes_dir = os.path.realpath(genomes_dir)
    start_run_time = time()
    # first create the empty index of k-mers
    index = BuildIndex(k)

    processed_dirs = 0
    total_dirs = len(os.listdir(genomes_dir))
    for species_dir in os.listdir(genomes_dir):
        time_start = time()
        species_dir = os.path.normpath(os.path.join(genomes_dir, species_dir))
        Clusterize(species_dir, k, G, F, J, L, index, max_seeds)
        processed_dirs += 1
        
        # When finished, set index again to empty
        # arrays in order to be used by next pangenome
        index = index.fromkeys(index, None)
        
        runtime = time() - start_run_time
        time_end = time() - time_start
        print "Pangenome calculated in {} seconds".format(time_end)
        print "Processed {} of {} species in {} seconds".format(processed_dirs, total_dirs, runtime)
        
        index["start_offset"] = 0

def main():
    import sys

    args = parse_args()
    genome_dir = args.genome_dir
    # scanned = args.scanned
    G = int(args.G)
    k = int(args.k)
    F = int(args.F)
    if not F > k:
        F = k + 1
        sys.stderr.write("f must be greater than k. Adjusting f to {}\n".format(F))
    J = int(args.J)
    L = int(args.L)
    N = float(args.N)
    max_seeds = int(args.max_seeds)
    recursive = args.recursive
    out_orphan = args.out_orphan
    
    if recursive:
        ProcessGenomesDir(genome_dir, k, G, F, J, L, N, max_seeds)
    elif not recursive:
        ProcessDir(genome_dir, k, G, F, J, L, N, max_seeds, out_orphan)
    else:
        assert False

if __name__ == "__main__":
    main()



