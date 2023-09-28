"""Ex2.py - BLAST (online)"""

import argparse
from Bio.Blast import NCBIWWW
from Bio import SeqIO

def perform_blast(sequence, output_file):
    """Realiza un BLAST de una secuencia y guarda los resultados en un archivo."""
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence, format_type="Text")
    with open(output_file, "w", encoding="utf8") as output_handle:
        output_handle.write(result_handle.read())
    result_handle.close()

def main(input_fasta_file, output_file):
    """Aplica BLAST sobre secuencias de aminoácidos en formato FASTA y guarda los resultados en un único archivo."""
    input_sequences = list(SeqIO.parse(input_fasta_file, "fasta"))
    for sequence in input_sequences:
        perform_blast(sequence.seq, output_file)
    print(f"{len(input_sequences)} secuencia{'s' if len(input_sequences) > 1  else ''} en {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply BLAST to FASTA amino acid sequences and save the output.")
    parser.add_argument("--input", required=True, help="The input FASTA file name.")
    parser.add_argument("--output", required=False, help="The output file name.")
    args = parser.parse_args()
    INPUT_FAS_FILE = f"output/fastaV2_{args.input}.fas"
    if args.output is None:
        OUTPUT_FILE = f"output/blast_{args.input}.out"
    else:
        OUTPUT_FILE = f"output/{args.output}.out"
    main(INPUT_FAS_FILE, OUTPUT_FILE)
