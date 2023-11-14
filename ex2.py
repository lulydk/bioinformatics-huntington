"""Ex2.py - BLAST (online)"""

import argparse
from Bio.Blast import NCBIWWW
from Bio import SeqIO

def perform_blast(sequence, output_file, format_type="Text"):
    """Realiza un BLAST de una secuencia y guarda los resultados en un archivo."""
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence, format_type=format_type)
    with open(output_file, "w", encoding="utf8") as output_handle:
        output_handle.write(result_handle.read())
    result_handle.close()

def main(input_fasta_file, output_file, format_type="Text"):
    """Aplica BLAST sobre secuencias de aminoácidos en formato FASTA y guarda los resultados en un único archivo."""
    input_sequences = list(SeqIO.parse(input_fasta_file, "fasta"))
    for sequence in input_sequences:
        perform_blast(sequence.seq, output_file, format_type=format_type)
    print(f"{len(input_sequences)} secuencia{'s' if len(input_sequences) > 1  else ''} en {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply BLAST to FASTA amino acid sequences and save the output.")
    parser.add_argument("--input", required=True, help="The input FASTA file name.")
    parser.add_argument("--output", required=False, help="The output file name.")
    parser.add_argument("--format", required=False, help="The output file format (Text or XML).")
    args = parser.parse_args()
    INPUT_FAS_FILE = f"output/fastaV2_{args.input}.fas"
    if args.output is None:
        OUTPUT_FILE = f"output/blast_{args.input}.out"
    else:
        OUTPUT_FILE = f"output/{args.output}.out"
    OUTPUT_FORMAT = args.format if args.format is not None else "Text"

    main(INPUT_FAS_FILE, OUTPUT_FILE, OUTPUT_FORMAT)
