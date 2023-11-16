import argparse
import random
import primer3


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

def extract_primers_info(primers_dict):
    extracted_primers = []

    for i, left_primer_info in enumerate(primers_dict['PRIMER_LEFT']):
        right_primer_info = primers_dict['PRIMER_RIGHT'][i]

        # Verificar condiciones adicionales si es necesario
        if (
            left_primer_info['GC_PERCENT'] >= 50.0 and
            left_primer_info['GC_PERCENT'] <= 60.0 and
            right_primer_info['GC_PERCENT'] >= 50.0 and
            right_primer_info['GC_PERCENT'] <= 60.0 and
            left_primer_info['TM'] <= 67.0 and
            right_primer_info['TM'] <= 67.0 and
            left_primer_info['SEQUENCE'][0] not in {'c', 'g'} and
            right_primer_info['SEQUENCE'][-1] not in {'g', 'c'}
        ):
            extracted_primers.append((left_primer_info['SEQUENCE'], right_primer_info['SEQUENCE']))

    return extracted_primers
def design_specific_primers(sequence, num_primers=5):
    primers_found = []
    # Configuración de los parámetros para Primer3
    while len(primers_found) < num_primers:
        seq_args = {
                'SEQUENCE_ID': f'my_sequence_{len(primers_found) + 1}',
                'SEQUENCE_TEMPLATE': sequence,
                'PRIMER_OPT_SIZE': random.randint(18, 24),
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 24,
                'PRIMER_OPT_TM': random.uniform(55.0, 65.0),
                'PRIMER_MIN_TM': 55.0,
                'PRIMER_MAX_TM': 65.0,
                'PRIMER_MIN_GC': 50.0,
                'PRIMER_MAX_GC': 60.0,
                'PRIMER_PRODUCT_SIZE_RANGE': [(100, 300)],
        }
        # Diseño de primers con Primer3
        potential_primers = primer3.bindings.design_primers(seq_args, {})
        primers_info= extract_primers_info(potential_primers)
        primers_found.extend(primers_info)
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
    parser.add_argument("--output", required=False, help="The output txt file name for primers.")
    args = parser.parse_args()

    INPUT_GB_FILE = f"input/{args.input}.gb"
    # OUTPUT_PRIMERS_FILE_RAND = args.output or f"output/primers_output_rand{args.input}.txt"
    OUTPUT_PRIMERS_FILE_PRIMER3 = args.output or f"output/{args.output}.txt"

    my_sequence = get_sequence(INPUT_GB_FILE)
    # primers = generar_primers(my_sequence, num_primers=5)
    primers_3 = design_specific_primers(my_sequence,5)
    write_primers(primers_3, OUTPUT_PRIMERS_FILE_PRIMER3)
    print(f"Primers generados fueron escritos en '{OUTPUT_PRIMERS_FILE_PRIMER3}'.")
