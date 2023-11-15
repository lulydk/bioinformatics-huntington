import argparse
import random


def generar_primers(sequence, num_primers=5, longitud_min=18, longitud_max=24, gc_min=50, gc_max=60, tm_max=67):
    primers_found = []

    while len(primers_found) < num_primers:
        longitud = random.randint(longitud_min, longitud_max)
        start_index = random.randint(0, len(sequence) - longitud)

        primer_candidate = sequence[start_index:start_index + longitud]

        # Calcular %GC
        gc_content = (primer_candidate.count('g') + primer_candidate.count('c')) / float(longitud) * 100

        # Calcular temperatura de melting (Tm) con corrección de sal
        tm = 41 * (gc_content - 16.4) / float(longitud)

        # Verificar condiciones
        if gc_min <= gc_content <= gc_max and tm <= tm_max:
            # Evitar GC en extremos terminales
            if primer_candidate[0] not in {'c', 'g'} and primer_candidate[-1] not in {'g', 'c'}:
                primers_found.append(primer_candidate)

    return primers_found


def get_sequence(input_file):
    sequence = ''
    with open(input_file, "r") as archivo:
        start_reading = False
        for linea in archivo:
            if "ORIGIN" in linea:
                start_reading = True
                continue
            if start_reading:
                # Filtra solo las letras 'g', 'c', 't', 'a' y elimina otros caracteres
                letras_validas = ''.join(caracter for caracter in linea if caracter.lower() in {'g', 'c', 't', 'a'})
                sequence = ''.join([sequence, letras_validas])
    return sequence


def write_primers(primers, output_file):
    with open(output_file, 'w') as archivo:
        for i, primer in enumerate(primers, 1):
            archivo.write(f"{primer}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate primers from a DNA sequence.")
    parser.add_argument("--input", required=True, help="The input GenBank file name.")
    parser.add_argument("--output", required=False, help="The output file name for primers.")
    args = parser.parse_args()

    INPUT_GB_FILE = f"input/{args.input}.gb"
    OUTPUT_PRIMERS_FILE = args.output or f"output/primers_output_{args.input}.txt"

    my_sequence = get_sequence(INPUT_GB_FILE)
    primers = generar_primers(my_sequence, num_primers=5)

    write_primers(primers, OUTPUT_PRIMERS_FILE)
    print(f"Primers generados fueron escritos en '{OUTPUT_PRIMERS_FILE}'.")
