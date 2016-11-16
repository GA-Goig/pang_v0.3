#! /usr/bin/python
def RecordsSizeParse(infile):

    # Skip any text before the first record (e.g. blank lines, comments)
    with open(infile) as handle:
        while True:
            line = handle.readline()
            if line == "":
                return  # Premature end of file, or just empty?
            if line[0] == ">":
                seek_value = handle.tell() - len(line)
                record_lines = 1
                break

        while True:
            if line[0] != ">":
                raise ValueError(
                    "Records in Fasta files should start with '>' character")
            #title = line[0:].rstrip()
            lines = []
            line = handle.readline()
            while True:
                if not line:
                    break
                if line[0] == ">":
                    break
                lines.append(line.rstrip())
                record_lines += 1
                line = handle.readline()

            # Remove trailing whitespace, and any internal spaces
            # (and any embedded \r which are possible in mangled files
            # when not opened in universal read lines mode)

            seq = "".join(lines).replace(" ", "").replace("\r", "")
            yield seek_value, record_lines, len(seq)
            seek_value = handle.tell() - len(line)
            record_lines = 1

            if not line:
                return  # StopIteration


def SortMultiFasta(infile, outfile):

    records = []
    for seek, lines, length in RecordsSizeParse(infile):
        records.append( (seek, lines, length) )

    # Sort records by decreasing length
    records = sorted(records, key=lambda x:x[2])

    with open(outfile, "w") as outfh:
        with open(infile) as infh:
            while records:
                next_record = records.pop()
                seek_pos = next_record[0]
                record_lines = next_record[1]
                infh.seek(seek_pos)
                for i in range(record_lines):
                    line = infh.readline()
                    outfh.write(line)

def main():
    from sys import argv
    SortMultiFasta(argv[1], argv[2])

if __name__ == "__main__":
    main()