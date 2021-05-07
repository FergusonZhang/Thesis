# This program parses VCF file with FSATA file.
import argparse
import vcf


def file_parser(filename):
    filenames = filename.split(".")
    file_type = filenames[1]
    if file_type == "fas":
        print("This is a FASTA file.")
        sequence = []
        with open(filename) as f:
            current_sequence = ""
            for line in f:
                if line[0] != '>':
                    current_sequence = current_sequence + line.strip()
                else:
                    if current_sequence:
                        sequence.append(current_sequence)
                        current_sequence = ""
            if current_sequence:
                sequence.append(current_sequence)
        return sequence
    elif file_type == "vcf":
        print("This is a VCF file.")
        # TODO: handle the VCF file.
    else:
        print("Unrecognized type of input!", type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genome reader.')
    parser.add_argument(dest='file_name', help='Please input the genome filename.')
    args = parser.parse_args()
    file_parser(args.file_name)
