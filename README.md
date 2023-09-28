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

El resultado estará en output/output_NM_000311-5.fas