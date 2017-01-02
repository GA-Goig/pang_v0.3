
def AlignRecords(fasta, index, index_map, k, G, F, J, L, N, max_seeds,
    cluster_title, curr_cluster, curr_orphan, out_orphan, 
    clusters={}, orphans={}, cmapping={}, omapping={}):
    '''This function iterates over each record in the fasta file and performs
    the SeedAndExtend function over each one, returning coordinates of "core"
    alignments for each record and indexed sequence in a dict where record
    title is the key'''
    import sys
    from parse_utils import FastaParser, GetID, GetRecordGroups
    from coordinates_utils import SortCoordinates, JoinFragments
    from coordinates_utils import  MapAlignments, GetNewCoreSeqs
    from index_utils import ReindexRecord
    from seq_utils import NucleotideFreq
    from aligner import SeedAndExtend

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
                            cluster = str(curr_cluster)
                            # Format header to CORE_TITLE format
                            header = cluster_title + cluster
                            # Add new sequences to index and update index_map
                            index = ReindexRecord(header, k, index, index_map, new_seq)
                            # Update pangenome dictionary with new_core_seq
                            clusters[header] = new_seq
                            # Update mapping_dict with new groups
                            ref_coords = (1, len(new_seq), ref, "plus", 
                                         new_seq_start + 1, new_seq_end + 1)
                            cmapping[cluster] = [ref_coords]
                            curr_cluster += 1
                        elif out_orphan:
                            orphan = str(curr_orphan)
                            header = cluster_title + orphan
                            orphans[header] = new_seq
                            orphan_coords = (1, len(new_seq), ref, "plus", 
                                             new_seq_start + 1, new_seq_end + 1)
                            omapping[orphan] = [orphan_coords]
                            curr_orphan += 1

                else:
                    # If no alignment is produced, scanned_core_coords is evaluated
                    # as False since it contains an empty list of coordinates
                    # in that case all sequence is new
                    if len(seq) > L and NucleotideFreq(seq, "N") < N:
                        cluster =  str(curr_cluster)
                        header = cluster_title + cluster
                        index = ReindexRecord(header, k, index, index_map, seq)
                        # Update pangenome with a new_core_seq
                        clusters[header] = seq
                        # Update mapping dict for that new_seq
                        ref_coords = (1, len(seq), ref, "plus", 1, len(seq))
                        cmapping[cluster] = [ref_coords]
                        curr_cluster += 1
                        # If all sequence is new, there is no need to look which records
                        # each alignment maps to, so continue with following iteration of
                        # for loop
                        continue
                    elif out_orphan:
                        orphan = str(curr_orphan)
                        header = cluster_title + orphan
                        orphans[header] = seq
                        orphan_coords = (1, len(seq), ref, "plus", 1, len(seq))
                        omapping[orphan] = [orphan_coords]
                        curr_orphan += 1
                        continue
                    else:
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
                        cmapping[cluster].append(
                            (c_start, c_end, ref, strand, scan_start, scan_end ))

 
    if title == "FILE_EMPTY": # If no title, seq returned by parser title
                              # remains <<FILE_EMPTY>>
        sys.exit("Error: scanned fasta file seems to be empty")

    to_return = index, index_map, curr_cluster, curr_orphan, cluster_title, \
    clusters, orphans, cmapping, omapping 
    
    return to_return

def Clusterize(fasta, k, G, F, J, L, N, max_seeds, update, out_orphan):
    '''Initialize necessary files and variables and call main function 
    << AlignRecords >>
    '''
    from init_pang import RealPathFiles, BaseName
    from index_utils import BuildIndex, StoreIDX
    from pang_IO import WriteClusters, WriteMapping
    from time import time

    base_name = BaseName(fasta)
    files = RealPathFiles(fasta)
    cmapfile, clustfile, omapfile, orphanfile, idx, fastfile = files

    # If we are not updating build a new index and new index_map 
    # otherwise load them from .idx and .imap files
    if not update:
        index = BuildIndex(k)
        index_map = [ [], [] ]
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
    AlignRecords(fastfile, index, index_map, k, G, F, J, L, N, max_seeds, 
        cluster_title, curr_cluster, curr_orphan, out_orphan)
    # Write clusters and mappings that are stored in memory
    clusters, orphans, cmapping, omapping = idx_info[-4:]
    WriteClusters(clusters, clustfile)
    WriteMapping(cmapping, cmapfile)
    if out_orphan:
        WriteClusters(orphans, orphanfile)
        WriteMapping(omapping, omapfile)

    time_end = time() - time_start
    print "{} clustered in {} seconds".format(fasta, time_end)

    # print "Saving .idx file..."
    # StoreIDX(idx_info, idx)
    # print "Done!"