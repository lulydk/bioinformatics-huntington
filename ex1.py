"""Ex1.py - PROCESAMIENTO DE SECUENCIAS"""

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def translate_sequence(dna_sequence):
    """Función para traducir una secuencia de nucleótidos a una secuencia de aminoácidos."""
    return dna_sequence.translate()

def main(input_genbank_file, output_fasta_file):
    """Traduce features CDS de un archivo GenBank a una secuencia de aminoácidos en formato FASTA."""
    translated_sequences = []
    with open(input_genbank_file, "r", encoding="utf8") as genbank_handle:
        for record in SeqIO.parse(genbank_handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    nucleotide_sequence = feature.location.extract(record).seq      # obtengo secuencia de nucleótidos indicada por el record
                    amino_acid_sequence = translate_sequence(nucleotide_sequence)
                    fasta_record = SeqRecord(amino_acid_sequence, id=feature.qualifiers.get("gene")[0])     # creo un registro FASTA con la secuencia traducida
                    translated_sequences.append(fasta_record)
    with open(output_fasta_file, "w", encoding="utf8") as fasta_handle:
        SeqIO.write(translated_sequences, fasta_handle, "fasta")    # escribo la(s) secuencia(s) traducida(s) en un archivo FASTA
    print(f"{len(translated_sequences)} secuencia{'s' if len(translated_sequences) > 1  else ''} en {output_fasta_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate CDS features from GenBank file to FASTA amino acid sequences.")
    parser.add_argument("--input", required=True, help="The input GenBank file name.")
    parser.add_argument("--output", required=False, help="The output FASTA file name.")
    args = parser.parse_args()
    INPUT_GB_FILE = f"input/{args.input}.gb"
    if args.output is None:
        OUTPUT_FAS_FILE = f"output/output_{args.input}.fas"
    else:
        OUTPUT_FAS_FILE = f"output/{args.output}.fas"
    main(INPUT_GB_FILE, OUTPUT_FAS_FILE)
