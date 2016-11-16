Old Functions

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

def ListFasta(genome_dir_path, compressed=False):
    '''This function gets all files within genome_dir and keeps only those with
    fasta, fna or fa extension, returning a list with them sorted descending 
    size in bytes

    This way it is possible to take the biggest genome file in order to use it
    as first reference-indexed genome to begin genome alignments.

    It is intended to work with directories that hold several genome files, from
    2 or 3 to thousands of them, plus some necessary files like .mapping or 
    .cinfo
    '''
    import re
    import os
    from os.path import normpath, join
    
    # List all files within genome_dir
    files = os.listdir(genome_dir_path)
    # Compile regexp to match only fasta|fna|fa extension files
    if compressed:
        extension = re.compile("\.gz$")
    else:
        extension = re.compile("\.(fasta|fna|fa)$")
        # Keep only files with that extension
    fasta_files = filter(extension.search, files)
    
    fasta_files = [ normpath(join(genome_dir_path, fasta)) for fasta in fasta_files ]
    # Sort fasta_files in-place by ascending size so the last element will be 
    # the larger genome
    fasta_files.sort(key=os.path.getsize)

    return fasta_files

def LoadInfo(core_info):
    '''This function loads genus, species and taxid info
    about the species being clustered so proper file names and headers
    can be used'''

    with open(core_info) as infile:
        # Read genus, species and taxid info
        genus = infile.readline().rstrip()
        genus = genus.split("=")[1]
        species = infile.readline().rstrip()
        species = species.split("=")[1]
        taxid = infile.readline().rstrip()
        taxid = taxid.split("=")[1]

        # According to this info, build a header to use as common
        # part for cluster sequences
        cluster_title = ">taxid|" + taxid + "|@" + genus + "_" + species + "@"

    
    return (genus, species, cluster_title)

def WritePangenome(pgnome_dict, pgnome_file):
     # Open with append, to add new sequences
     del pgnome_dict["CURRENT"]
     del pgnome_dict["TITLE"]
     with open(pgnome_file, "w") as outfile:
        for record in pgnome_dict:
            new_seq = pgnome_dict[record]
            outfile.write(record + "\n")
            # The ruler is to write sequences
            # formated with a maximum of 100
            # columns
            ruler = 0
            for i in xrange(len(new_seq)):
                if ruler < 101:
                    outfile.write(new_seq[i])
                    ruler += 1
                else:
                    outfile.write("\n" + new_seq[i])
                    ruler = 0
            outfile.write("\n")

def WriteMapping(mapping_dict, mapping_file):
    try:
        with open(mapping_file, "w") as outfile:
            for cluster in mapping_dict:
                alignments = mapping_dict[cluster]
                for alignment in alignments:
                    c_start, c_end, acc, strand, seq_start, seq_end = alignment
                    score = len(alignments)
                    string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        cluster, c_start, c_end, acc, score, strand, seq_start, seq_end
                        )
                    outfile.write(string)
    except TypeError:
        print "Exception : {} : {}".format(cluster, mapping_dict[cluster])