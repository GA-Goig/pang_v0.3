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
    parser.add_argument("--work-on-disk", dest="wok", action="store_true",
        help="If enabled, Mondo will use much less RAM but will run significantly slower too")

    args = parser.parse_args()

    return args




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
    
    if args.wok:
        from pang.work_on_disk import Clusterize
    else:
        from pang.work_on_ram import Clusterize

    Clusterize(fasta, k, G, F, J, L, N, max_seeds, update, out_orphan)
    
if __name__ == "__main__":
    main()



