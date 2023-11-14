
from Bio.Blast import NCBIXML
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from Bio import SeqIO
import argparse

def extract_top_sequences(xml_file, output_file='top_10_sequences.fasta'):
    # Parse the BLAST report file
    result_handle = open(xml_file)
    blast_records = NCBIXML.parse(result_handle)

    # Extract the top 10 sequences
    top_sequences = []
    for record in blast_records:
        for alignment in record.alignments[:10]:
            sequence_id = alignment.title.split()[1]  # Extracting sequence ID from the title
            sequence_length = alignment.length
            sequence = alignment.hsps[0].sbjct

            top_sequences.append((sequence_id, sequence_length, sequence))

    # Write top 10 sequences to a FASTA file
    with open(output_file, 'w') as fasta_file:
        for index, (sequence_id, sequence_length, sequence) in enumerate(top_sequences, start=1):
            fasta_file.write(f">Sequence_{index} | ID: {sequence_id} | Length: {sequence_length}\n")
            fasta_file.write(f"{sequence}\n")

    print(f"Top 10 sequences exported to {output_file}")


def perform_msa_with_pairwise2(fasta_file, output_file='msa_result.txt'):
    # Parse the FASTA file
    sequences = list(SeqIO.parse(open(fasta_file, 'r'), 'fasta'))

    # Use the first sequence as the reference
    original_seq = sequences[0].seq

    # Perform pairwise alignment for each sequence with the original sequence
    alignments = []
    for sequence_record in sequences[1:]:
        sequence = sequence_record.seq
        alignment = pairwise2.align.globalxx(original_seq, sequence)
        alignments.append(format_alignment(*alignment[-1]))

    # Write alignments to the output file
    with open(output_file, 'w') as save_file:
        save_file.write('\n'.join(alignments))

    print(f"Multiple sequence alignment completed. Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract Top 10 from Blast Report and perform Multiple Sequence Alignment using pairwise2')
    #add argument for blast xml report
    parser.add_argument('-x', metavar='XML_FILE', help='Input BLAST XML report (default = blast_report.xml)', default='output/blast_report.xml')
    parser.add_argument('-i', metavar='FASTA_FILE', help='Generated FASTA file (default = top_10_sequences.fasta)', default='output/top_10_sequences.fasta')
    parser.add_argument('-o', metavar='MSA_FILE', help='Output MultiSequence Alignment file (TXT) (default = msa_results.txt)', default='output/msa_results.txt')
    args = parser.parse_args()

    xml_file_path = args.x
    fasta_file_path = args.i
    output_msa_file = args.o

    extract_top_sequences(xml_file_path, fasta_file_path)

    perform_msa_with_pairwise2(fasta_file_path, output_msa_file)
