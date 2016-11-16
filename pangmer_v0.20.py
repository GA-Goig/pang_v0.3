#!/usr/bin/pypy
# -*- coding: utf-8 -*-


def parse_args():
    '''Parse arguments given to script'''

    import argparse

    parser = argparse.ArgumentParser(description="Given two genomes, get both "\
        "common and specific regions to both genomes following a kmer-based algorithm")

    parser.add_argument("-in", dest="fasta", required=True, metavar="Input sorted multi fasta file")
    
    parser.add_argument("-g", metavar="jumping Gap between k-mers. DEFAULT: 6", dest="G", 
        default=6, help="Gap between k-mers. From inf to -k+1. (e.g G = 0 if no"\
        " gaps are desired. G = -K+1 if an overlap of 1 base between k-mers is desired")
    
    parser.add_argument("-k", metavar="k-mer length. DEFAULT: 12", dest="k", default=12, 
        help="k-mer length. Avoid use k > 12 in desktop computers")
    
    parser.add_argument("-f", metavar="min alignment length. DEFAULT: 100", dest="F", default=100,
        help="Distance checked by k-mer jumping from a seed k-mer to be considered as" \
        " a significant alignment")
    
    parser.add_argument("-j", metavar="max distance to combine fragments. DEFAULT: 10", 
        dest="J", default=10, help="Maximum distance to join significant alignments "\
        "from k-mer jumping")
    
    parser.add_argument("--max-seeds", metavar="maximum seeds per kmer. DEFAULT: 20", 
                        dest="max_seeds", default=20, help="Maximum seeds of a k-mer"\
                        " in order to seed alignments. This limit allows for a dramatic "\
                        "speedup in repetitive regions")
    
    parser.add_argument("-l, --length", dest="L", metavar="min cluster length. DEFAULT: 500", 
        default=500, help="Minimum length of a sequence in order to be considered as a new Cluster")
    
    parser.add_argument("-n", dest="N", metavar="max percentage of N per cluster. DEFAULT: 0.3",
     default=0.3, help="Maximum percentage of ambiguous nucleotides in clusters. If a new"\
     " cluster exceeds this maximum, it is discarded")

    parser.add_argument("-o", dest="out_orphan", action="store_true",
        help="Output orphan clusters that did not pass L and N filters"\
        "(.orphan and .omapping)")
    parser.add_argument("-u", dest="update", action="store_true", 
        help="Use -u if new sequences have to be added to previous clustered files")

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

def AlignRecords(files, index, index_map, k, G, F, J, L, N, max_seeds, 
    cluster_title, curr_cluster, curr_orphan, out_orphan):
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
    from pang.init_pang import BaseName
    from pang.pang_IO import UpdateMapping, UpdateClusters

    # Get clustering files names from input fasta name
    cmapping, clusters, omapping, orphans, idx, fasta = files

    # NOTE that each time a coordinate is going to be appended to a mapping dict
    # we sum +1 so we switch from programming to human coordinates
    title = "FILE_EMPTY"
    with open(fasta) as handle:
        for title, seq in FastaParser(handle):
            # Get the ref from current record
            ref = GetID(title)
            if seq:
                # First produce alignments with indexed sequences and retrieve
                # coordinates of these alignments
                alignment_coordinates = SeedAndExtend(seq, index, k, G, F, max_seeds)
                if title == "M02918:57:000000000-APPEG:1:2118:19555:24274":
                    print "alignment_coordinates = {}".format(alignment_coordinates)
                if alignment_coordinates:
                    # Sort alignment coordinates
                    sorted_coordinates = SortCoordinates(alignment_coordinates)
                    # And Join fragments that are as close as parameter J
                    joined_coords = JoinFragments(sorted_coordinates, J)
                    # Get new sequences that do not produce core alignments
                    # This new sequences will be added to a new GROUP so first
                    # update the global variable with the new group
                    if title == "M02918:57:000000000-APPEG:1:2118:19555:24274":
                        print "joined_coords = {}".format(joined_coords)
                    new_seqs = GetNewCoreSeqs(joined_coords, seq, L, N)
                    for new_seq, new_seq_start, new_seq_end, orphan in new_seqs:
                        if title == "M02918:57:000000000-APPEG:1:2118:19555:24274":
                            print "new_seq = {}".format(new_seq)
                            print "new_seq_start = {}".format(new_seq_start)
                            print "new_seq_end = {}".format(new_seq_end)
                        if not orphan:
                            cluster = "Cluster_" + str(curr_cluster)
                            # Format header to CORE_TITLE format
                            header = cluster_title + cluster
                            # Add new sequences to index and update index_map
                            index = ReindexRecord(header, k, index, index_map, new_seq)
                            # Update .clusters  file
                            UpdateClusters(clusters, header, new_seq)
                            # Update .cmapping file
                            ref_coords = (1, len(new_seq), ref, "plus", 
                                         new_seq_start + 1, new_seq_end + 1)
                            UpdateMapping(cmapping, cluster, ref_coords)
                            curr_cluster += 1
                        elif out_orphan:
                            orphan = "Orphan_" + str(curr_orphan)
                            header = cluster_title + orphan
                            # Update .orphan file
                            UpdateClusters(orphans, header, new_seq)
                            # Update .omapping file
                            orphan_coords = (1, len(new_seq), ref, "plus", 
                                             new_seq_start + 1, new_seq_end + 1)
                            UpdateMapping(omapping, orphans, orphan_coords)
                            curr_orphan += 1

                else:
                    # If no alignment is produced, scanned_core_coords is evaluated
                    # as False since it contains an empty list of coordinates
                    # in that case all sequence is new
                    if len(seq) >= L and NucleotideFreq(seq, "N") < N:
                        cluster = "Cluster_" + str(curr_cluster)
                        header = cluster_title + cluster
                        index = ReindexRecord(header, k, index, index_map, seq)
                        # Update .cluster file
                        UpdateClusters(clusters, header, seq)
                        # Update .cmapping 
                        ref_coords = (1, len(seq), ref, "plus", 1, len(seq))
                        UpdateMapping(cmapping, cluster, ref_coords)
                        curr_cluster += 1
                        # If all sequence is new, there is no need to look which records
                        # each alignment maps to, so continue with following iteration of
                        # for loop
                        continue
                    elif out_orphan:
                        orphan = "Orphan_" + str(curr_orphan)
                        header = cluster_title + orphan
                        # Update .orphans
                        UpdateClusters(orphans, header, seq)
                        # Update .omapping
                        orphan_coords = (1, len(seq), ref, "plus", 1, len(seq))
                        UpdateMapping(omapping, orphan, orphan_coords)
                        curr_orphan += 1
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
                        # Update .cmapping
                        coords = (c_start, c_end, ref, strand, scan_start, scan_end)    
                        UpdateMapping(cmapping, cluster, coords)

    if title == "FILE_EMPTY": # If no title, seq returned by parser title
                              # remains <<FILE_EMPTY>>
        sys.exit("Error: scanned fasta file seems to be empty")

    return index, index_map, curr_cluster, curr_orphan, cluster_title

def Clusterize(fasta, k, G, F, J, L, N, max_seeds, update, out_orphan):
    '''Initialize necessary files and variables and call main function 
    << AlignRecords >>
    '''
    from pang.init_pang import RealPathFiles, BaseName
    from pang.index_utils import BuildIndex, StoreIDX
    from time import time

    base_name = BaseName(fasta)
    files = RealPathFiles(fasta)
    # Get idx and imap file names to store/load serialized objects
    idx = files[-2]

    # If we are not updating build a new index and new index_map 
    # otherwise load them from .idx and .imap files
    if not update:
        index = BuildIndex(k)
        index_map = [ [] ,[] ]
        curr_cluster = 1
        curr_orphan = 1
        cluster_title = "{}@".format(base_name)
        # Truncate possible files from previous runs
        for file in files[:-1]: # Last file is the input fasta
            fh = open(file, "w")
            fh.close()
    else:       
        index, index_map, curr_cluster, curr_orphan, cluster_title = LoadIDX(idx)

    time_start = time()
    print "Clustering {}...".format(fasta)
    idx_info = \
    AlignRecords(files, index, index_map, k, G, F, J, L, N, max_seeds, 
        cluster_title, curr_cluster, curr_orphan, out_orphan)
    time_end = time() - time_start
    print "{} clustered in {} seconds".format(fasta, time_end)

    # print "Saving .idx file..."
    # StoreIDX(idx_info, idx)
    # print "Done!"

def main():
    import sys

    args = parse_args()
    fasta = args.fasta
    # scanned = args.scanned
    G = int(args.G)
    k = int(args.k)
    F = int(args.F)
    if not F > k:
        F = k + 1
        sys.stderr.write("f must be greater than k. f set to {}\n".format(F))
    J = int(args.J)
    L = int(args.L)
    N = float(args.N)
    max_seeds = int(args.max_seeds)
    update = args.update
    out_orphan = args.out_orphan
    
    Clusterize(fasta, k, G, F, J, L, N, max_seeds, update, out_orphan)

if __name__ == "__main__":
    main()



