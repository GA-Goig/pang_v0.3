def UpdateClusters(clusters, header, new_seq):
    '''Write a new cluster to clusters_file'''

    # Open with append mode
    with open(clusters, "a") as outfh:
        WriteSeq(outfh, header, new_seq)

def UpdateMapping(mapping, cluster, alignment):
    '''Updates a mapping file with new records for a cluster alignment or new
    clusters'''
    c_start, c_end, acc, strand, seq_start, seq_end = alignment
    string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                                        cluster, 
                                                        c_start, 
                                                        c_end, 
                                                        acc, 
                                                        strand, 
                                                        seq_start, 
                                                        seq_end
                                                        )
    with open(mapping, "a") as outfh:
        outfh.write(string)

        
def WriteNewCoreSeqs(new_core_seq, core_file, header):
    '''Writes new sequences added to core genome in the core genome file
  
    This function takes some GLOBAL variables:

    current_group is an int value of the new group to be generated to map te gi
    identifier where new sequences come from

    CORE_TITLE is a string of the common part of the header of all sequences
    comprising the core genome
    
    For example, building a Mycobacterium tuberculosis core genome

    >taxid|1773| Mycobacterium_tuberculosis Group:1

    Group 1 will correspond to some gi in an additional file called mapping_file
    So in this example, if Group 1 is mapping to gi's X, Y, Z, you can tell this
    sequence is a Mycobacterium tuberculosis sequence and, more specifically, 
    common to sequences with gene identifiers X, Y, Z
    
    So the common part, that is stored <<CORE_TITLE>> would be:

    >taxid|1773| Mycobacterium tuberculosis Group:

    '''
    
    with open(core_file, "a") as handle:
        WriteSeq(handle, header, new_core_seq)

def WriteSeq(handle, header, seq, ruler=101):
    '''This function takes a fasta header and sequence and writes
    sequence in fasta format with ruler as a value for formating
    lines'''

    handle.write(">" + header + "\n")
    rul = 0
    for i in xrange(len(seq)):
        if rul < ruler:
            handle.write(seq[i])
            rul += 1
        else:
            handle.write("\n" + seq[i])
            rul = 0
    handle.write("\n")