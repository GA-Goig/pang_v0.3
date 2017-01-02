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
        if new_seq:
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
    print "start, end_coord = {},  {}".format(start_coord, end_coord)
    print "map_start, end   = {},  {}".format(map_start, map_end)
    print "record = {}".format(record)
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

def JoinFragments(sorted_coords, J):
    '''This function takes a list of sorted tuples of tuples of
    (start,end) coordinates associated to coordinates where the
    scanned sequence has aligned in index and joins them if a
    distance less than J separate them
    
    In the following example for J=20

    e.g 

    INPUT = [
             ((0,100), [(100,200), (600, 500), (1000000, 1000100)]),
             ((120,500), [(220,600), (480, 100)])
            ]
    
    OUTPUT = [ 
             ((0,100) [(1000000,1000100)])
             ((0,500),[(100,600), (600,100)])
             ]
    '''
    from copy import copy


    scoords = copy(sorted_coords)
    # Join fragments that are no more distant that J value
    i = 0
    # Singletons is a list that will hold coordinates for fragments that have
    # not been joined so in a last step they are added to joined fragments
    # and both, joined and unjoined fragments are returned by the function
    singletons = []
    while i < len(scoords) - 1:
        # Take contiguous t uples of coordinates
        tuple_A = scoords[i][0]
        tuple_B = scoords[i+1][0]
        index_coords_A = scoords[i][1]
        index_coords_B = scoords[i+1][1]
        joined = False

        # Join vectors keep track of what coordinates have been joined, 0 means
        # that coordinate has not been joined while 1 means joining with another
        # pair of coordinates
        join_vector_A = [0]*len(index_coords_A)
        join_vector_B = [0]*len(index_coords_B)

        new_index_coords = []
        # Iterate over A index coordinates trying to join them with B index coords
        for j in range(len(index_coords_A)):
            index_A = index_coords_A[j]
            # Check if are coords of a reverse alignment, strand will hold True
            # if reverse and False if forward
            strand_A = AreReverseCoords(index_A)
            for k in range(len(index_coords_B)):
                index_B = index_coords_B[k]
                strand_B = AreReverseCoords(index_B)
                # Only try joining if they are in same strand
                if strand_A == strand_B:
                    # Calc distance according to its orientation 
                    if strand_A: # If they are reverse coordinates
                        distance = index_A[1] - index_B[0]
                    else: # If they are forward
                        distance = index_B[0] - index_A[1]
                    # If distance is less or equal than user defined J then join
                    if distance <= J and distance >= 0:
                        # Update join vectors
                        join_vector_A[j] = 1
                        join_vector_B[k] = 1
                        joined = True
                        # Get new joined coords
                        new_index_start = index_A[0]
                        new_index_end = index_B[1]
                        # And update new index coords list
                        new_index_coords.append( (new_index_start, new_index_end) )
        
        if not joined: # If no join has been produced, pass to next tuple
            i += 1
        else: # If any join has been produced, update coords
            new_start = tuple_A[0]
            new_end = tuple_B[1]
            scoords[i] = ( (new_start, new_end), new_index_coords)
            del scoords[i+1]
            # Update singletons list with fragments that have not been joined
            no_join_A = []
            for j in range(len(index_coords_A)):
                if join_vector_A[j] == 0:
                    no_join_A.append(index_coords_A[j])
            if no_join_A:
                singletons.append((tuple_A, no_join_A))

            no_join_B = []
            for k in range(len(index_coords_B)):
                if join_vector_B[k] == 0:
                    no_join_B.append(index_coords_B[k])
            if no_join_B:
                singletons.append((tuple_B, no_join_B))

    # merge with singletons
    merged = scoords + singletons
    # resort and return
    return sorted(merged, key=lambda x:x[0])

def AreReverseCoords(coords):
    '''Return True if coords are pointing out a reverse alignment, otherwise
    return False

    e.g (100, 200) --> False
        (580, 320) --> True

    '''
    if coords[0] > coords[1]:
        return True
    else:
        return False

    # This line should not be reached as a pair of alignment coordinates should
    # never be equal
    assert False

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