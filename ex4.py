import subprocess

# Paso 1: Análisis de ORFs
def analisis_orfs(input_file, output_file):
    subprocess.run(["getorf", "-sequence", input_file, "-outseq", output_file])

# Paso 2: Descarga de motivos de PROSITE
def descargar_prosite_motivos(propositedir="prosite"):
    subprocess.run(["wget","https://ftp.expasy.org/databases/prosite/prosite.dat"])
    subprocess.run(["wget","https://ftp.expasy.org/databases/prosite/prosite.doc"])
    subprocess.run(["chmod","775"])
    subprocess.run(["mkdir", propositedir])
    subprocess.run(["mv", "prosite.dat", propositedir])
    subprocess.run(["mv", "prosite.doc", propositedir])
    subprocess.run(["prosextract", propositedir])

# Paso 3: Análisis de dominios de secuencias de aminoácidos
def analisis_dominios(input_file, output_file):
    subprocess.run(["patmatmotifs", "-sequence", input_file, "-outfile", output_file])

if __name__ == "__main__":

    # Ejecutar el análisis de ORFs
    analisis_orfs("./input/NM_001388492.1.gb", "./output/secuencia_orfs.fasta")

    # Descargar motivos de PROSITE
    descargar_prosite_motivos()

    # Ejecutar el análisis de dominios de secuencias de aminoácidos
    analisis_dominios("./output/fastaV2_NM_001388492.1.fas", "output/analis_dominios.txt")
