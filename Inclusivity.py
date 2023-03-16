"""
coding=utf-8
author = Alejandro González

"""
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from filtros import process
from ProbeOk import classifier


def compare_seq(query_seq, subject_seq):
    replace = []
    aux = 0
    seq = ""

    for x, y in zip(query_seq, subject_seq):
        if len(query_seq) == len(subject_seq):
            if x == y:
                replace.append("-")
            else:
                replace.append(y)
                aux += 1
        seq = "".join(replace)
    return seq, aux


def inclusivity(filename, oligos, excel_name):
    df = pd.read_csv(filename, sep="\t", header=None)  # Archivo .txt, resultado del blast
    df.columns = ['Oligo', 'Referencia', 'Rstart', 'Rend', 'Ostart', 'Oend', 'Mismatch', 'Kingdom', 'Group', 'Group2',
                  'Strand', 'Evalue', 'QuerySeq', 'SubjectSeq']
    df['maxL'] = df.Oend - df.Ostart + 1  # Calculamos el tamaño para comprobar si hay acortamientos
    # sentido = pd.DataFrame([df['Oligo'], df['Referencia'], df['Strand']]).T

    # Llamamos a los primers
    records = list(SeqIO.parse(oligos, "fasta"))

    olidf = pd.DataFrame({'LENGTH': [len(records[0].seq), len(records[1].seq), len(records[2].seq)],
                          'Oligo': ['FORWARD', 'REVERSE', 'PROBE']})

    datamerged = df.merge(olidf, on='Oligo')  # Unimos los datos que hemos obtenido del blast con los primers

    datamerged = datamerged.drop_duplicates(subset=['Oligo', 'Referencia'])
    datamerged = datamerged.reset_index()

    #  Esta col cuenta los mismatch del principio o final mirando si coinciden longitudes
    datamerged['Acortamiento'] = datamerged.LENGTH - datamerged.maxL

    newdata3, finalls, nF, nR, nP, newdatafyr3, finallsfyr, nF2, nR2, selectedAlignments_df = process(datamerged)

    nP.columns = ['Referencia', 'Pstart', 'P-Acortamiento', 'P-Mismatch']
    nR.columns = ['Referencia', 'Rstart', 'R-Acortamiento', 'R-Mismatch']

    selectedAlignments_df.columns = ['Referencia', 'Fstart', 'F-Acortamiento', 'F-Mismatch', 'Rstart', 'Pstart']

    # Dividimos las referencias en los diferentes primers (Forward, Reverse y Probe)
    nFF = newdata3[newdata3['Oligo'] == 'FORWARD']
    nRR = newdata3[newdata3['Oligo'] == 'REVERSE']
    nPP = newdata3[newdata3['Oligo'] == 'PROBE']

    nFF = nFF.rename(columns={'Rstart': 'Fstart'})
    nPP = nPP.rename(columns={'Rstart': 'Pstart'})

    # Esta tabla corresponde a la columa aling en la pestaña raw, el objetivo es comprobar la homologia de los primers
    # con el resultado obtenido despues de hacer el blast

    """
    Recorremos a la vez nuestro primer con el tipo de primer (Forward, Reverse o Probe).
    El segundo bucle, compara nuestro primer con el resultado obtenido al realizar el blast para cada referencia. 
    Esta comparacion da como resultado la homología de nuestro primer con la secuencia que hemos comprado.

    Por ejemplo:
        1.- CGATGAGGCTATTCCGACTAGGT	
            CGATGAGGCTATTCCGACTAGGT 
            -----------------------
            Como se puede ver, al coincidir nuestro primer completamente con la secuencia que hemos comparado el 
            resultado es 100% homologo.

        2.- ATGAGGCTATTCCGACTAGGT	
            ATGAGCATATTCCGACTAGGT
            -----CA--------------

            En este segundo caso, hay dos nn que no coinciden, por lo que el primer no es 100% homologo, sino que 
            tiene 2 mm.
    """

    # Recorrer los todos los primers y con el Oligo
    counter = 0
    total = []
    tamano = []

    for i, j, b in zip(datamerged['QuerySeq'], datamerged['Oligo'], tqdm(range(len(datamerged)))):
        z = 0
        seq, aux = compare_seq(i, datamerged['SubjectSeq'][counter])
        total.append(seq)
        counter += 1

        if j == 'FORWARD':
            z = (aux / len(records[0].seq)) * 100
            z = 100 - z
        elif j == 'REVERSE':
            z = (aux / len(records[1].seq)) * 100
            z = 100 - z
        elif j == 'PROBE':
            z = (aux / len(records[2].seq)) * 100
            z = 100 - z

        tamano.append(z)

    # Se añaden las dos columnas obtenidas (Alineamiento y porcentaje de homologia) a la tabla total
    datamerged['Aling'] = total
    datamerged['Accuracy'] = tamano

    # Evitar alineamientos multiples
    datamerged = datamerged.drop_duplicates(subset=['Oligo', 'Referencia'])

    # Esta tabla corresponde a la pestaña data
    final = selectedAlignments_df.merge(nP[['Referencia', 'Pstart', 'P-Acortamiento', 'P-Mismatch']], on='Referencia')
    final = final.merge(nR[['Referencia', 'Rstart', 'R-Acortamiento', 'R-Mismatch']], on='Referencia')
    final = final.drop("Rstart_x", axis=1)  # Eliminar columnas innecesarias
    final = final.drop("Pstart_x", axis=1)  # Eliminar columnas innecesarias
    final.rename(columns={"Rstart_x": "Rstart", "Pstart_x": "Pstart"}, inplace=True)

    final = final.merge(nFF[['Referencia', 'Fstart', 'Group', 'Kingdom', 'Group2']], on=['Referencia'])
    final = final.merge(nPP[['Referencia', 'Pstart', 'Group', 'Kingdom', 'Group2']], on=['Referencia'])
    final = final.merge(nRR[['Referencia', 'Rstart', 'Group', 'Kingdom', 'Group2']], on=['Referencia'])
    final = final.apply(classifier, axis='columns')

    final['resta'] = final.Rstart - final.Fstart_x  # Resta para determinar si hay acortamiento
    final = final.fillna('Indeterminado')
    final = final.reindex(columns=['Referencia', 'Fstart_x', 'F-Acortamiento', 'F-Mismatch', 'Pstart', 'P-Acortamiento',
                                   'P-Mismatch', 'Rstart', 'R-Acortamiento', 'R-Mismatch', 'Group', 'Kingdom', 'Group2',
                                   'resta', 'ProbeOk'])
    final = final.drop_duplicates()

    # Esta tabla corresponde a la pestaña data_filter
    final_filter = final

    # Aplicar filtros para quitar aquellas secuencias que no nos interesan
    # Tamaño del amplicon de 500 pb maximo y los primers en orden correcto
    final_filter = final_filter.loc[(final_filter.ProbeOk == 'Yes') & (final_filter.resta < 500) &
                                    (final_filter.resta > -500)]

    final_filter = final_filter.reindex(columns=['Referencia', 'Fstart_x', 'F-Acortamiento', 'F-Mismatch', 'Pstart',
                                                 'P-Acortamiento', 'P-Mismatch', 'Rstart', 'R-Acortamiento',
                                                 'R-Mismatch', 'Group', 'Kingdom', 'Group2', 'resta', 'ProbeOk'])

    # print("Va a descargar las secuencias fasta")

    # Eliminar duplicados para evitar alineamientos multiples
    final_filter = final_filter.drop_duplicates()

    # Corresponde a la pestaña referencias_repetidas
    referencias_repetidas = []
    localizador = []
    counter = []

    # En el caso de que haya alineamientos multiples, calcula el número de veces que se repite dicho alineamiento
    for i, a, w in zip(final_filter.Referencia, final_filter.Group2, tqdm(range(len(final_filter.Referencia)))):
        contador = 0

        for j in final_filter.Referencia:

            if i == j:
                contador += 1

        if contador > 1:
            referencias_repetidas.append(i)
            localizador.append(a)
            counter.append(contador)

    referencias_repetidas = pd.DataFrame(referencias_repetidas, columns=['Referencia'])
    referencias_repetidas['Nombre'] = localizador
    referencias_repetidas['Num_veces'] = counter
    referencias_repetidas = referencias_repetidas.drop_duplicates()

    # Esta tabla corresponde a la pestaña Accuracy
    Accuracy = datamerged

    # Separar la lista por primers
    f = list(Accuracy['Referencia'].loc[(Accuracy.Oligo == 'FORWARD')])
    r = list(Accuracy['Referencia'].loc[(Accuracy.Oligo == 'REVERSE')])
    p = list(Accuracy['Referencia'].loc[(Accuracy.Oligo == 'PROBE')])

    # Tener solo en cuenta las referencias que se encuentran en las tres listas, ya que puede haber referencias a
    # secuencias que solo se encuentran en uno de los tres primers. Ese caso no amplificaria, por lo que no se tiene
    # en cuenta
    t = [i for i in f if i in r and i in p]

    Acc = Accuracy

    # Eliminar aquellas referencias que no cumplen la condicion anterior
    index = [Acc.index[j] for i, j in zip(Acc['Referencia'], range(len(Acc))) if i not in t]
    Acc = Acc.drop(index)

    # Eliminar las columnas que no necesitamos en esta tabla
    Ax = Acc.drop(columns=['Rstart', 'Rend', 'Ostart', 'Oend', 'Mismatch', 'Strand', 'Evalue', 'QuerySeq', 'SubjectSeq',
                           'maxL', 'LENGTH', 'Acortamiento', 'Accuracy'])
    Ax = Ax.sort_values(by='Referencia')

    # Calcula el porcentaje para cada una de las referencias que cumplen la condicion
    new_colum = []

    for ref, c in zip(Ax['Referencia'], tqdm(range(len(Ax)))):
        suma = 0
        tam = 0

        ref_rows = Ax[Ax['Referencia'] == ref]

        for _, row in ref_rows.iterrows():
            suma += row['Aling'].count('-')
            lens = row['Oligo']

            if lens == "FORWARD":
                tam += len(records[0].seq)
            elif lens == "REVERSE":
                tam += len(records[1].seq)
            elif lens == "PROBE":
                tam += len(records[2].seq)

        total = (suma / tam) * 100
        new_colum.append(total)

    Ax['Total_Accuracy'] = new_colum  # Añade el porcentaje a la tabla

    Ax_final = Ax.drop_duplicates(subset=['Referencia'])
    Ax_final = Ax_final.drop(columns=['Oligo', 'Aling'])

    # Esta tabla corresponde a la pestaña resume
    final2 = final.drop_duplicates(subset=['Referencia'])
    final2 = final2.groupby(['Kingdom', 'Group', 'Group2']).Referencia.count()

    # Esta tabla corresponde a la pestaña ProbeOk
    final3 = final.drop_duplicates(subset=['Referencia'])
    final3 = final3[(final3.ProbeOk == 'Yes')].groupby(['Kingdom', 'Group', 'Group2']).Referencia.count()

    # Esta tabla corresponde a la pestaña mm
    final3mm_filter = final
    # Aplicamos los filtros del tamaño del amplicon y tener los primers en el orden correcto
    final3mm_filter = final3mm_filter.loc[(final3mm_filter.ProbeOk == 'Yes') &
                                          (final3mm_filter.resta < 500) &
                                          (final3mm_filter.resta > -500)]
    # Contamos el numero de veces que se repite cada alineamiento
    final3mm = final3mm_filter.groupby(['Kingdom', 'Group', 'Group2', 'F-Acortamiento', 'F-Mismatch', 'P-Acortamiento',
                                        'P-Mismatch', 'R-Acortamiento', 'R-Mismatch', 'resta']).Referencia.count()

    # Esta tabla corresponde a la pestaña fyR
    # Para las tablas que siguen, solo tiene en cuenta el Forward y el Reverse
    finallsfyr = finallsfyr.drop_duplicates(subset=['Referencia'])
    finalfyr = finallsfyr.merge(nF2[['Referencia', 'Rstart', 'Acortamiento']], on='Referencia')
    finalfyr = finalfyr.merge(nR2[['Referencia', 'Rstart', 'Acortamiento']], on='Referencia')
    finalfyr.columns = ['Referencia', 'Fstart', 'F-MM', 'Rstart', 'R-MM']
    finalfyr = finalfyr.merge(newdatafyr3[['Referencia', 'Group', 'Kingdom', 'Group2']], on='Referencia')
    finalfyr = finalfyr.drop_duplicates()
    finalfyr['resta'] = finalfyr.Rstart - finalfyr.Fstart
    finalfyr = finalfyr.fillna('Indeterminado')

    # Mismo protocolo que las tablas anteriores, pero solo para Forward y Reverse
    # Esta tabla corresponde a la pestaña fyR_filtred
    finalfyr_filtred = finalfyr[(finalfyr.resta < 500) & (finalfyr.resta > -500)]  # No tiene en cuenta el ProbeOk
    # porque no hay Probe

    # Esta tabla corresponde a la pestaña fyR_group
    finalfyr2 = finalfyr_filtred.groupby(['Kingdom', 'Group', 'Group2']).Referencia.count()

    # Esta tabla corresponde a la pestaña fyR_mm
    finalfyrmm = finalfyr.groupby(['Kingdom', 'Group', 'Group2', 'F-MM', 'R-MM']).Referencia.count()

    print("Creando el excel ...")

    # Crear Excel
    datamerged.to_excel(excel_name, sheet_name='raw')
    with pd.ExcelWriter(excel_name, mode='a', engine="openpyxl") as writer:
        final.to_excel(writer, sheet_name='data')
        final_filter.to_excel(writer, sheet_name='data_filter')
        referencias_repetidas.to_excel(writer, sheet_name='referencias_repetidas')
        Ax_final.to_excel(writer, sheet_name='Accuracy')
        final2.to_excel(writer, sheet_name='resume')
        final3.to_excel(writer, sheet_name='ProbeOk')
        final3mm.to_excel(writer, sheet_name='mm')
        finalfyr.to_excel(writer, sheet_name='fyR')
        finalfyr_filtred.to_excel(writer, sheet_name='fyR_filtred')
        finalfyr2.to_excel(writer, sheet_name='fyR_group')
        finalfyrmm.to_excel(writer, sheet_name='fyR_mm')

    print("Ha acabado")
