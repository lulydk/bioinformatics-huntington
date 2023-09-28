"""Ex1.py - PROCESAMIENTO DE SECUENCIAS"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def translate_sequence(dna_sequence):
    """Función para traducir una secuencia de nucleótidos a una secuencia de aminoácidos."""
    return dna_sequence.translate()

INPUT_GB_FILE = "input/NM_000311-5.gb"
OUTPUT_FAS_FILE = "output/output_NM_000311-5.fas"

translated_sequences = []

with open(INPUT_GB_FILE, "r", encoding="utf8") as genbank_handle:
    for record in SeqIO.parse(genbank_handle, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                nucleotide_sequence = feature.location.extract(record).seq      # obtengo secuencia de nucleótidos indicada por el record
                amino_acid_sequence = translate_sequence(nucleotide_sequence)
                fasta_record = SeqRecord(amino_acid_sequence, id=feature.qualifiers.get("gene")[0])     # creo un registro FASTA con la secuencia traducida
                translated_sequences.append(fasta_record)

with open(OUTPUT_FAS_FILE, "w", encoding="utf8") as fasta_handle:
    SeqIO.write(translated_sequences, fasta_handle, "fasta")    # escribo la(s) secuencia(s) traducida(s) en un archivo FASTA

print(f"{len(translated_sequences)} secuencia{'s' if len(translated_sequences) > 1  else ''} en {OUTPUT_FAS_FILE}")
