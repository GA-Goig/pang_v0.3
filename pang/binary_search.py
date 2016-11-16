# This module takes a multi fasta file and outputs
# its records sorted by decreasing length in another
# filehander

def SortMultiFasta(infile, outfile):
    from bisect import insort_right
    from pang.pang_IO import WriteSeq
    from pang.parse_utils import FastaParser

    records_size = []
    with open(infile) as infh:
        for title, seq in FastaParser(infh):
            records_size.append( (title, len(seq)) )

    # Sort the list
    records_size = sorted(records_size, key=lambda x:x[1])

    while records_size:
        bigest_record = records_size.pop()[0]
        with open(infile) as infh:
            with open(outfile, "a") as outfh:
                for title, seq in FastaParser(infh):
                    if title == bigest_record:
                        WriteSeq(outfh, title, seq)
                        bigest_record = records_size.pop()[0]

def main():
    from sys import argv
    SortMultiFasta(argv[1], argv[2])

if __name__ == "__main__":
    main()