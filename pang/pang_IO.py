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

def WriteMapping(mapping, mapfile):
    '''Dump clusters info in clusters file'''
    with open(mapfile, "w") as infh:
        for cluster in mapping:
            alignments = mapping[cluster]
            for alignment in alignments:
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
                infh.write(string)

def WriteClusters(clusters, clustfile):
    '''Dump mapping info in mapping file'''
    with open(clustfile, "w") as outfh:
        for seqid in clusters:
            seq = clusters[seqid]
            WriteSeq(outfh, seqid, seq)

