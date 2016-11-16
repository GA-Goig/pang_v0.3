def BaseName(fasta):
    '''Try to get a base name from fasta file name provided as input

    e.g. In = "Mycobacterium_tuberculosis.sorted.fasta"
        Out = "Mycobacterium_tuberculosis.sorted"
    '''
    import re

    search = re.search("(.*)(.fasta|.fas|.fna)", fasta)
    if search:
        # Get rid off fasta extension
        base_name_path = search.groups()[0]
        # Get rid off path
        base_name = base_name_path.split("/")
        return base_name[-1]
    else:
        return fasta

def BaseFiles(fasta):
    ''' Return base_files names from input fasta file
    base_name.cmapping
    base_name.clusters
    base_name.omapping
    base_name.orphan
    .base_name.idx
    '''

    base_name = BaseName(fasta)

    cmapping = "{}.cmapping".format(base_name)
    clusters = "{}.clusters".format(base_name)
    omapping = "{}.omapping".format(base_name)
    orphan = "{}.orphan".format(base_name)
    idx = ".{}.idx".format(base_name)

    return cmapping, clusters, omapping, orphan, idx

def RealPathFiles(fasta):
    '''Get full path files
    e.g /home/Jon_Doe/M_tuberculosis.fasta instead of M_tuberculosis.fasta
    '''
    from os.path import realpath
    
    fasta_file = realpath(fasta)

    # Get path for same directory where input fasta is stored
    file_items = fasta_file.split("/")
    working_path = "/".join(file_items[:-1])

    files = BaseFiles(fasta)
    real_files = ["{}/{}".format(working_path, fname) for fname in files]

    real_files.append(fasta_file)

    return real_files


