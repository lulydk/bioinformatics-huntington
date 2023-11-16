"""Ex1.py - PROCESAMIENTO DE SECUENCIAS"""

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_longest_proteins_by_orf(record: SeqRecord):
    """Retorna un diccionario con las distintas proteínas obtenidas al variar el Open Reading Frame (ORF)"""
    orf_proteins = {}
    strand_sequences = [(+1, record.seq), (-1, record.seq.reverse_complement())]
    for strand, nucleotides in strand_sequences:
        for frame in range(1, 4):
            frame_offset = (4 - frame) % 3 # number of nucleotides to skip from the beggining of the sequence
            record_to_translate = nucleotides[frame_offset:] # extract subsequence
            record_to_translate = record_to_translate[:3 * (len(record_to_translate) // 3)] # seq length as mult of 3 (required by biopython)
            longest_protein_length = 0
            for protein in record_to_translate.translate().split('*'): # * signals the end (stop codon)
                protein = protein[protein.find('M'):] # M signals the start
                if protein.startswith('M') and len(protein) > longest_protein_length:
                    orf = strand * frame # orf in [1, 2, 3, -1, -2, -3]
                    orf_proteins[f'+{orf}' if orf > 0 else str(orf)] = protein
                    longest_protein_length = len(protein)
    return orf_proteins

def find_possible_orf(orf_proteins: dict):
    """Encuentra el ORF más probable seleccionando la proteína más larga."""
    most_probable_orf = 0
    longest_protein_length = 0
    for orf, protein_sequence in orf_proteins.items():
        protein_length = len(protein_sequence)
        if protein_length > longest_protein_length:
            longest_protein_length = protein_length
            most_probable_orf = orf
    return most_probable_orf


def main(input_genbank_file, output_fasta_file):
    """Obtiene una secuencia de aminoácidos en formato FASTA a partir de un archivo GenBank."""
    with open(input_genbank_file, "r", encoding="utf8") as genbank_handle:
        for record in SeqIO.parse(genbank_handle, "genbank"):
            proteins = get_longest_proteins_by_orf(record)
            best_orf = find_possible_orf(proteins)
            with open(output_fasta_file, "w", encoding="utf8") as fasta_handle:
                fasta_record = SeqRecord(proteins[best_orf], id=record.id, description=f"Using ORF {best_orf}")
                SeqIO.write(fasta_record, fasta_handle, "fasta")    # escribo la(s) secuencia(s) traducida(s) en un archivo FASTA

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate from GenBank file to FASTA amino acid sequences.")
    parser.add_argument("--input", required=True, help="The input GenBank file name.")
    parser.add_argument("--output", required=False, help="The output FASTA file name.")
    args = parser.parse_args()
    INPUT_GB_FILE = f"input/{args.input}.gb"
    if args.output is None:
        OUTPUT_FAS_FILE = f"output/fastaV2_{args.input}.fas"
    else:
        OUTPUT_FAS_FILE = f"output/{args.output}.fas"
    main(INPUT_GB_FILE, OUTPUT_FAS_FILE)
