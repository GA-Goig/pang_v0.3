def SortCoordinates(coordinates, reverse=False):
    '''This function simply sort a list of pairs os coordinates in tuples
    by ascending value of the start coordinate (first number of each tuple)

    If reverse=True, then SortCoordinates expect reverse_complement aligned
    coords to sort


    e.g  INPUT  = [(500, 100), (1000, 450), (80, 50)]
         OUTPUT = [(1000, 450), (500, 100), (80, 50)]



    '''
    if reverse:
        return sorted(coordinates, key=lambda x:x[1], reverse=True)
    else:
        return sorted(coordinates, key=lambda x:x[0])

def GetNewCoreSeqs(joined_coords, seq, L, N):
    '''GENERATOR: Given a list of sorted/fragment_joined coordinates of core 
    alignments, and the sequence that produced those alginments, Yield sequences
    that DO NOT BELONG TO THE CORE, (opposite seqs to coordinates provided), so 
    they could be reindexed and added as new core sequences
    
    Coordinates of regions aligned marked with ( )  

              FIRST                                      REMAINING
    Input =  |------(-------)--------(---)-------(-----)--------|
    Output = |-----|         |------|     |-----|       |-------| 
    
    Yield sequence that have at least length L
    '''
    from pang.seq_utils import NucleotideFreq
    
    # joined_coords are tuples pairs of scanned_coordinates with an
    # asociated list of indexed coordinates they aligned with, get
    # only scanned coordinates in order to get new sequences
    scanned_coords = [coords[0] for coords in joined_coords]
    # Check first if there are new seqs in the beginning  of the seq
    first_coords = scanned_coords[0] # Take first coordinates
    first_start = first_coords[0]
    # Boolean to check if new_seq pass L and N filters
    orphan = False
    # Check if first aligned coord is at the beginning of scanned seq
    if first_start != 0:
        # If first aligned coord is not at the beginning of seq
        # take first part of sequence as new_seq (as in example above)
        new_seq = seq[0 : first_start]
        if len(new_seq) < L or NucleotideFreq(new_seq, "N") > N:
            orphan = True
        new_seq_start = 0
        new_seq_end = first_start
        yield (new_seq, new_seq_start, new_seq_end, orphan)

    for i in xrange(len(scanned_coords) - 1):
        orphan = False
        tuple_A = scanned_coords[i] # First pair of coordinates (previous)
        tuple_B = scanned_coords[i+1] # Next pair of coordinates
        previous_end = tuple_A[1]
        next_start = tuple_B[0]
        new_seq = seq[previous_end + 1 : next_start]
        if len(new_seq) < L or NucleotideFreq(new_seq, "N") > N:
            orphan = True
        new_seq_start = previous_end
        new_seq_end = next_start
        yield (new_seq, new_seq_start, new_seq_end, orphan)  

    # Take the remaining sequence after last pair of coordinate if it didn't produce
    # any alignment
    # Check if the end of sequence is within a region aligned
    last_coords = scanned_coords[-1]
    end_coord = last_coords[1]
    orphan = False
    if end_coord < len(seq) - 1: # Then the end of seq is not within an alignment
        new_seq = seq[end_coord + 1 : len(seq)]
        if len(new_seq) < L or NucleotideFreq(new_seq, "N") > N:
            orphan = True
        new_seq_start = end_coord
        new_seq_end = len(seq)
        yield (new_seq, new_seq_start, new_seq_end, orphan)  

def MapCoordinates(index_map, coords, reverse=False):
    '''This function takes the index_map and a pair of start end coordinates
    and returns all records that coordinates belong to

    Index = |--------|    |-------|    |-------|    |-------|
            0   R1   10  20  R2   30  40   R3  50  50   R4  60
  
                |--------------------------|
                5                          45

    Given start_coord = 5, end_coord = 45, result will be 
    [(R1,5,10), (R2, 0, 10), (R3,0,5)]

    If reverse_complement coordinates are provided, then start and 
    end coords are fliped out
    '''

    records = []
    if not reverse:
        start_coord = coords[0]
        end_coord = coords[1]
    else:
        start_coord = coords[1]
        end_coord = coords[0]

    # Iterate over each pair of coordinates of index_map
    map_records = index_map[0]
    map_coords = index_map[1]
    for i in range(len(map_coords)):
        map_start = map_coords[i][0] # Start index coordinate for record i
        map_end = map_coords[i][1] # End index coordinate for record i
        if start_coord < map_end:
            record = map_records[i]
            if start_coord >= map_start:
                start = start_coord - map_start
            elif start_coord < map_start:
                start = 0
            if end_coord > map_end:
                end = map_end - map_start
            elif end_coord <= map_end:
                end = end_coord - map_start
                if not reverse:
                    records.append((record, start, end))
                else:
                    records.append((record, end, start))
                return records
            if not reverse:
                records.append((record, start, end))
            else:
                records.append((record, end, start))
    # This line should not be reached. Only would be if end_coord is > than map_end
    # and that should not be possible
    assert False

def MapAlignments(joined_coords, index_map):
    ''' This function takes a list of tuples of coordinates describing alignments
    each tuple contains both, a tuple of coordinates of the scanned sequence that
    produced alignments in the index, and a list of tuples of coordinates 
    describing which index regions have been aligned with. This function returns
    which clusters are comprised by that coordinates and which cluster_coordinates
    are implied in the alignment. EXAMPLE:

    Input  = [ ( (100, 200), [(500,600), (100200, 100300)] )]
    Output = [ ( (100, 200), [("Cluster_1, 0, 100"), ("Cluster_51, 555, 655")] )]
    '''

    records_map = []
    for scan_coords, index_coords in joined_coords:
        records_matched = []
        # Split in forward and reverse index alignments coordinates
        Fw_coords, Rv_coords = SplitFwRv(index_coords)
        # Map coordinates for Fw_coords
        for coords in Fw_coords:
            records = (MapCoordinates(index_map, coords))
            for record in records:
                records_matched.append(record)
        # And map coordinates for Rv_coords
        for coords in Rv_coords:
            records = (MapCoordinates(index_map, coords, reverse=True))
            for record in records:
                records_matched.append(record)

        records_map.append( (scan_coords, records_matched) )

    return records_map

def JoinCoordinates(sorted_coords, J, reverse=False):
    '''This function takes a list of tuples of start, end coordinates sorted
    by ascending value and try yo join them if end coordinate of one tuple is
    less than J nucleotides of distance that start coordinate of the next tuple

    eg.  IN  = [ (100, 200), (225, 333), (1000, 2000)]; J=50
         OUT = [ (100, 333), (1000, 2000)]
    '''
    from copy import copy

    scoords = copy(sorted_coords)
    # Join fragments that are no more distant that J value
    i = 0
    while i < len(scoords) - 1:
        # Take contiguous t uples of coordinates
        tuple_A = scoords[i]
        tuple_B = scoords[i+1]

        # Take coordinates of each tuple
        start_coord_A = tuple_A[0]
        end_coord_A = tuple_A[1]
        start_coord_B = tuple_B[0]
        end_coord_B = tuple_B[1]
        distance = start_coord_B - end_coord_A # Then calculate distance
        # If sorted_coords describe a reverse_complement alignment, then 
        # distance will be always negative, so get the absolute value
        if reverse:
            distance = -distance
        if distance <= J: # If distance is lower than value J
            new_start = start_coord_A
            new_end = end_coord_B
            # Update second tuple with new coordinates and delete first tuple            
            scoords[i] = (new_start, new_end)  
            del scoords[i+1]  
        else:
            i += 1

    return scoords


def JoinFragments(sorted_coords, J):
    '''This function takes a list of sorted tuples of tuples of
    (start,end) coordinates associated to coordinates where the
    scanned sequence has aligned in index and joins them if a
    distance less than J separate them
    '''
    from copy import copy

    scoords = copy(sorted_coords)
    # Join fragments that are no more distant that J value
    i = 0
    while i < len(scoords) - 1:
        # Take contiguous t uples of coordinates
        tuple_A = scoords[i][0]
        tuple_B = scoords[i+1][0]
        index_coords_A = scoords[i][1]
        index_coords_B = scoords[i+1][1]

        # Take coordinates of each tuple
        start_coord_A = tuple_A[0]
        end_coord_A = tuple_A[1]
        start_coord_B = tuple_B[0]
        end_coord_B = tuple_B[1]
        distance = start_coord_B - end_coord_A # Then calculate distance         
        if distance <= J: # If distance is lower than value J
            new_start = start_coord_A
            new_end = end_coord_B
            new_index_coords = index_coords_A + index_coords_B
            # Split new index coordinates in two variables depending if they
            # describe alignments in forward or reverse strand
            forward_icoords, reverse_icoords = SplitFwRv(new_index_coords)
            # sort new forward and reverse icoords 
            forward_sorted = SortCoordinates(forward_icoords)
            reverse_sorted = SortCoordinates(reverse_icoords, reverse=True)
            # And join those separated less than J nucleotides
            forward_joined = JoinCoordinates(forward_sorted, J)
            reverse_joined = JoinCoordinates(reverse_sorted, J, reverse=True)
            # Now join forward and reverse coordinates again to return along
            # with scanned sequence
            joined_new_icoords = forward_joined + reverse_joined
            # Update second tuple with new coordinates and delete first tuple            
            scoords[i] = ( (new_start, new_end), joined_new_icoords)   
            del scoords[i+1]  
        else:
            i += 1

    return scoords

def SplitFwRv(coordinates):
    '''This function takes a list of tuples of coords and split it in
    two lists, one of forward alignments coordinates and the other one
    in reverse_complement alignments

    e.g  INPUT  = [(500,600), (1000,900)]
         OUTPUT =  ( Fw=[(500,600)], Rv=[(1000,900)] )
    '''

    Fw = []
    Rv = []
    for coords in coordinates:
        start = coords[0]
        end = coords[1]
        if start < end:
            Fw.append(coords)
        else:
            Rv.append(coords)

    return Fw, Rv
            
def CompactMap(mapping):
    '''This function takes a mapping dict and reduces its number of
    records if different cluster records actually come from same
    region aligned eg:
    
    INPUT:

    Cluster_2       1399    1536    NZ_BAYP01000052.1       5903    strand  1422    1992
    Cluster_2       1537    1674    NZ_BAYP01000052.1       5903    strand  1422    1992
    Cluster_2       1688    1969    NZ_BAYP01000052.1       5903    strand  1422    1992

    OUTPUT

    Cluster_2       1399    1969    NZ_BAYP01000052.1       5901    strand  1422    1992

    ACTUAL INPUT:

    mapping = { 
                "Cluster_2" : (1399, 1536, NZ_BAYP01000052.1, 5903, strand, 1422, 1992),
                ...
                ...
                }
    '''

    # Dict for final result initialized with same records and empty lists
    compact_mapping = mapping.fromkeys(mapping, [])
    # For each cluster
    for cluster in mapping:
        # Take all records
        records = mapping[cluster]
        # Begin from first record
        curr_record = records[0]
        # Store current cluster start, end, ref, sequence start, end for this record
        curr_cstart, curr_cend, curr_ref, cstrand, curr_sstart, curr_send = \
        curr_record
        # Now iterate over each record to see if they can be collapsed starting
        # from the second
        for next_record in records[1:]:
            next_cstart, next_cend, next_ref, nstrand, next_sstart, next_send = \
            next_record
            # If next record comes from same sequence accession (ref)
            if curr_ref == next_ref:
                # And from same sequence region
                if curr_sstart == next_sstart and curr_send == next_send:
                    # Both records in the same cluster are actually together, so
                    # collapse its coordinates in a single record
                    curr_cend = next_cend
                else:
                    # Although same ref, if they do not come from same region, they
                    # are actually separated records. So add current record to compact
                    record = (curr_cstart, curr_cend, curr_ref, cstrand, 
                              curr_sstart, curr_send)
                    compact_mapping[cluster].append(record)
                    # And current record is now the one we were processing
                    curr_cstart, curr_cend, curr_ref, cstrand, curr_sstart, curr_send = \
                    next_record
            else:
                # If they do not come from same ref, they are separated records,
                # so ad current record to compact
                record = (curr_cstart, curr_cend, curr_ref, cstrand, curr_sstart, 
                          curr_send)
                compact_mapping[cluster].append(record)
                # And current record is now the one we were processing
                curr_cstart, curr_cend, curr_ref, cstrand, curr_sstart, curr_send = \
                next_record
        # Append last record
        record = (curr_cstart, curr_cend, curr_ref, cstrand, curr_sstart, curr_send)
        compact_mapping[cluster].append(record)

    return compact_mapping
        