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

    jump = k + G
    gapped_kmer_gen = GappedKmerGenerator(sequence, k_start, jump, k)
    current_index_coord = seed_coordinate
    kmer = gapped_kmer_gen.next()

    # If alignment is being produced in reverse_complement, move "jump" backwards
    if reverse:
        jump = -jump
    while kmer in index:
        kmer_index_coords = index[kmer] # Take coords for next k+G gapped kmer

        if kmer_index_coords:
            current_index_coord += jump # Move k+G bases of index
            # When both coordinates don't coincide they are not contiguous
            if current_index_coord in kmer_index_coords:
                # If both coordinates coincide, get and check next gapped kmer
                kmer = gapped_kmer_gen.next()
            else:
                current_index_coord -= jump
                if not reverse:
                    alignment_length = (current_index_coord - seed_coordinate) + k 
                else:
                    alignment_length = (current_index_coord - seed_coordinate) - k

                return alignment_length
        else:
            current_index_coord -= jump
            if not reverse:
                alignment_length = (current_index_coord - seed_coordinate) + k 
            else:
                alignment_length = (current_index_coord - seed_coordinate) - k
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