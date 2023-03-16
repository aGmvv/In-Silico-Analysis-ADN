import time
import os
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio import Entrez


def LNA(filename, oligos, excel_name, patogeno):
    tic = time.perf_counter()
    df = pd.read_csv(filename, sep="\t", header=None)  # Archivo .txt, resultado del blast
    df.columns = ['Oligo', 'Referencia', 'Rstart', 'Rend', 'Ostart', 'Oend', 'Mismatch', 'Kingdom', 'Group', 'Group2',
                  'Strand', 'Evalue', 'QuerySeq', 'SubjectSeq']

    records = list(SeqIO.parse(oligos, "fasta"))

    # ------------------------------------------------------------------------------------------ #
    #    Este apartado descarga de NCBI las secuencias que cumplen con los filtros determinados  #
    #    anteriormente, y las guarda en un arhivo .fasta                                         #
    # ------------------------------------------------------------------------------------------ #

    Entrez.email = "agonzalez@certest.es"
    lista_id = df['Referencia']
    grupo = df['Group2']

    lista = pd.concat([lista_id, grupo], axis=1)
    lista = lista.drop_duplicates(subset=['Referencia', 'Group2'])

    interF = df.query('Oligo == "FORWARD"')
    dmF = interF.groupby(['Referencia']).Oligo.count()
    dmF = pd.DataFrame(dmF)
    dmF = dmF.reset_index()

    interR = df.query('Oligo == "REVERSE"')
    dmR = interR.groupby(['Referencia']).Oligo.count()
    dmR = pd.DataFrame(dmR)
    dmR = dmR.reset_index()

    regs = df.Referencia.unique()
    regs = pd.DataFrame(regs, columns=['Referencia'])
    newdata = regs.merge(dmF, on='Referencia', how='left')
    newdata = newdata.merge(dmR, on='Referencia', how='left')

    newdata.columns = ['Referencia', 'FORWARD', 'REVERSE']

    newdata2 = newdata[(newdata.FORWARD >= 1) & (newdata.REVERSE >= 1)]
    ls = list(newdata2.Referencia)

    finalls = pd.DataFrame(ls, columns=['Referencia'])
    newdata3 = df.query('Referencia in @ls')

    nF = newdata3[newdata3['Oligo'] == 'FORWARD']
    nR = newdata3[newdata3['Oligo'] == 'REVERSE']

    # nF = nF[['Referencia', 'Rstart', 'Mismatch']]
    # nR = nR[['Referencia', 'Rstart', 'Mismatch']]

    fasta = str(patogeno) + '.fasta'
    fw = open(fasta, 'w')

    # Recorre la lista de referencias y si corresponde al patogeno especificado, descarga
    # la secuencia de NCBI.
    for i, j, c in zip(lista.Referencia, lista.Group2, tqdm(range(len(lista_id)))):
        if patogeno in j or patogeno == j:
            hd1 = Entrez.efetch(db="nucleotide", id=[i], rettype='fasta')
            seq = SeqIO.read(hd1, 'fasta')
            SeqIO.write(seq, fw, 'fasta')

    # Cierra el archivo una vez acabado
    fw.close()
    os.getcwd()
    # ------------------------------------------------------------------------------------------ #

    print("FINISH")
    toc = time.perf_counter()
    print(f"Code in {toc - tic:0.4f} seconds")
