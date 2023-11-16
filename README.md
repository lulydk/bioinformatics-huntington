# TP Introducción a la Bioinformática

Para clonar el repositorio:

```bash
git clone https://github.com/lulydk/bioinformatics-huntington.git
cd bioinformatics-huntington
```

## Ejercicio 1

Traduce la región de codificación de un gen (CDS) de un archivo de GenBank a una secuencia de aminoácidos (usando el código genético estándar) según el reading frame especificado, y escribe la salida a un archivo en formato FASTA.

### Uso

Por línea de comando:

```bash
python3 ex1.py --input='NCBI_REF_SEQ' [ --output='OUTPUT_FILENAME' ]
```

De no especificar nombre del archivo de salida, se utiliza el mismo código provisto en el input.

### Ejemplo

Utilizando el [Homo sapiens prion protein (PRNP), transcript variant 1, mRNA - Nucleotide - NCBI (nih.gov)](https://www.ncbi.nlm.nih.gov/nuccore/NM_000311.5):

```bash
python3 ex1.py --input='NM_000311-5'
```

El resultado estará en `output/fastaV2_NM_000311-5.fas`

## Ejercicio 2

Aplica BLAST (online vía BioPython) a una secuencia de aminoácidos en formato FASTA, y escribe el resultado en un archivo de texto plano.

### Uso

Por línea de comando:

```bash
python3 ex2.py --input='NCBI_REF_SEQ' [ --output='OUTPUT_FILENAME' ]
```

De no especificar nombre del archivo de salida, se utiliza el mismo código provisto en el input.

### Ejemplo

Utilizando el archivo creado en el punto anterior:

```bash
python3 ex2.py --input "NM_000311-5"
```

El resultado estará en `output/blast_NM_000311-5.out`

## Ejercicio 3

Por línea de comando:

```bash
python3 ex3.py -x <XML_FILE> -i <FASTA_FILE> -o <MSA_FILE>
```

## Ejercicio 4

Por línea de comando:

```bash
python3 ex4.py --input=<ARCHIVO_FASTA> --output=<ANALISIS_DOM>
```

## Ejercicio 5

Por línea de comando:

```bash
python3 ex5.py --input=<ARCHIVO_GENBANK> --output=<PRIMERS>
```
