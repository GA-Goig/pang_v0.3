 
def KmerGenerator(sequence, k_start, k):
    '''GENERATOR: Given a sequence and a kmer length, yield k-mers from sequence
    moving one base at a time from k_start coordinate'''
    
    k_end = k_start + k
    while True:
        if k_end <= len(sequence):
            kmer = sequence[k_start : k_end]
            yield kmer
            k_start += 1
            k_end += 1
        else:
            yield 0

def GappedKmerGenerator(sequence, k_start, k, G):
    '''GENERATOR: Given a sequence, a kmer length and a gap value G, yield 
    k-mers from sequence moving k + G bases at a time from k_start coordinate'''
    
    jump = k+G
    k_start += jump
    k_end = k_start + k
    while True:
        if k_end <= len(sequence):
            kmer = sequence[k_start : k_end]
            yield kmer
            k_start += jump
            k_end += jump
        # If end of sequence has been reached, yield the offset between
        # end of sequence and last k_end position to check ends
        else:
            offset = k_end - len(sequence)
            yield offset


def SkipAmbiguous(sequence, k_start, non_ambiguous, min_non_ambiguous):
    '''This function is used when a k-mer with ambiguous nucleotides has been 
    found during SeedAndExtend. In that case, it takes the start position of
    that kmer in a sequence, and starts checking from that positions ambiguous
    nucleotides. When a number of CONSECUTIVE non-ambiguous nucleotides is
    checked, it returns position of last checked nucleotide. That position
    will be used as a new k_start parameter to call again SeedAndExtend.

    This "filter" has been implemented since simply calling SeedAndExtend
    recursively causes Python stop preventing an overflow in the name space

    non_ambiguous is the set { "A", "T", "G", "C" } passed by default to 
    SeedAndExtend'''

    count = 0
    seq_length = len(sequence)
    while count < min_non_ambiguous:
        if k_start < seq_length:
            nt = sequence[k_start]
            if nt in non_ambiguous:
                count += 1 # When that nucleotides is A, T, G or C
                k_start += 1  # Keep track how many nucleotides checked
            else:
                count = 0 # Restart counter when is an ambiguous nucleotide
                k_start += 1 # Keep track how many nucleotides checked
        # When while loop ends, min_non_ambiguous contiguous nucleotides have been
        # found. Then, we want SeedAndExtend start from first of that nucleotides
        # and that corresponds to current k_start - min_non_ambiguous
        else:
            # If max length of sequence has been reached, return last position
            # checked (which is k_start)
            return k_start
    
    k_start -= min_non_ambiguous

    return k_start

def NucleotideFreq(sequence, nucleotide):
    '''Given a sequence and a nucleotide, calculate the nucleotide frequence'''

    length = len(sequence)
    count = 0
    for nt in sequence:
        if nt == nucleotide:
            count += 1

    freq = count / float(length)
    
    return freq