import warnings
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


def aligmentselection(Fw_df, Rv_df, Pr_df1, Pr_df2, ls):
    new_df = pd.DataFrame(Fw_df)
    new_df = new_df.rename(columns={'Rstart': 'Fstart'})

    #  A partir de la lista con las referencias finales, se crean diccionarios de la siguiente manera:
    dR = {}  # dR['Referencia'] = (lista de posibles alineamientos → todas las posiciones de inicio que
    #  tiene el Rv en cada secuencia)
    for ref, c in zip(ls, tqdm(range(len(ls)))):
        dR[ref] = list(Rv_df['Rstart'].loc[Rv_df['Referencia'] == ref])

    dP1 = {}  # dP['Referencia'] = (lista de posibles alineamientos → todas las posiciones de inicio que tiene la Pr
    # en cada secuencia)
    dP2 = {}

    for ref, c in zip(ls, tqdm(range(len(ls)))):
        dP1[ref] = list(Pr_df1['Rstart'].loc[Pr_df1['Referencia'] == ref])
    for ref, c in zip(ls, tqdm(range(len(ls)))):
        dP2[ref] = list(Pr_df2['Rstart'].loc[Pr_df2['Referencia'] == ref])
    #  Comprobación de que el amplicon entre Fw y Rv tiene el tamaño adecuado
    new_df['Rstart'] = np.nan

    for ref, c in zip(ls, tqdm(range(len(ls)))):  # Se recorre la lista de referencias
        Fstarts = list(new_df['Fstart'].loc[new_df['Referencia'] == ref])

        for st in Fstarts:  # Para cada referencia, se recorre la si lista de posibles alineamientos del forward
            RvStart_temp = [Rstart for Rstart in dR[ref] if abs(Rstart - st) < 500]

            if len(RvStart_temp) >= 1:  # Si para un Fw hay varios Rv posibles, se elige el primero
                new_df['Rstart'].loc[new_df['Fstart'] == st] = RvStart_temp[0]  # Se añade la posicion start del
                #  Fw a una nueva col

    new_df = new_df.dropna()
    #  Comprobación de que la sonda está en medio entre Fw y Rv
    new_df['Pstart'] = np.nan

    for ref, c in zip(ls, tqdm(range(len(ls)))):
        Fstarts = list(new_df['Fstart'].loc[new_df['Referencia'] == ref])
        RVstarts = list(new_df['Rstart'].loc[new_df['Referencia'] == ref])
        tuple_FwRv = tuple(zip(Fstarts, RVstarts))

        for tupl in tuple_FwRv:
            fw = tupl[0]
            rv = tupl[1]
            PrStart_temp = [Rstart for Rstart in dP1[ref] or dP2[ref] if (fw < Rstart < rv) or (rv < Rstart < fw)]

            if len(PrStart_temp) >= 1:
                new_df['Pstart'].loc[new_df['Fstart'] == fw] = PrStart_temp[0]

    new_df = new_df.dropna()

    return new_df


def process(datamerged):
    #  Para Forward, Reverse y Sonda
    interF = datamerged.query('Oligo == "FORWARD"')
    dmF = interF.groupby(['Referencia']).Oligo.count()
    dmF = pd.DataFrame(dmF)
    dmF = dmF.reset_index()

    interR = datamerged.query('Oligo == "REVERSE"')
    dmR = interR.groupby(['Referencia']).Oligo.count()
    dmR = pd.DataFrame(dmR)
    dmR = dmR.reset_index()

    interP1 = datamerged.query('Oligo == "PROBE1"')
    dmP1 = interP1.groupby(['Referencia']).Oligo.count()
    dmP1 = pd.DataFrame(dmP1)
    dmP1 = dmP1.reset_index()

    interP2 = datamerged.query('Oligo == "PROBE1"')
    dmP2 = interP2.groupby(['Referencia']).Oligo.count()
    dmP2 = pd.DataFrame(dmP2)
    dmP2 = dmP2.reset_index()

    regs = datamerged.Referencia.unique()
    regs = pd.DataFrame(regs, columns=['Referencia'])
    newdata = regs.merge(dmF, on='Referencia', how='left')
    newdata = newdata.merge(dmR, on='Referencia', how='left')
    newdata = newdata.merge(dmP1, on='Referencia', how='left')
    newdata = newdata.merge(dmP2, on='Referencia', how='left')
    newdata.columns = ['Referencia', 'FORWARD', 'REVERSE', 'PROBE1', 'PROBE2']

    newdata2 = newdata[(newdata.FORWARD >= 1) & (newdata.REVERSE >= 1) & (newdata.PROBE1 >= 1) & (newdata.PROBE2 >= 1)]
    ls = list(newdata2.Referencia)

    finalls = pd.DataFrame(ls, columns=['Referencia'])
    newdata3 = datamerged.query('Referencia in @ls')

    nF = newdata3[newdata3['Oligo'] == 'FORWARD']
    nR = newdata3[newdata3['Oligo'] == 'REVERSE']
    nP1 = newdata3[newdata3['Oligo'] == 'PROBE1']
    nP2 = newdata3[newdata3['Oligo'] == 'PROBE2']

    nF = nF[['Referencia', 'Rstart', 'Acortamiento', 'Mismatch']]
    nR = nR[['Referencia', 'Rstart', 'Acortamiento', 'Mismatch']]
    nP1 = nP1[['Referencia', 'Rstart', 'Acortamiento', 'Mismatch']]
    nP2 = nP2[['Referencia', 'Rstart', 'Acortamiento', 'Mismatch']]

    #  Seleccionar alineamientos correctos comprobando
    selectedAlignments_df = aligmentselection(nF, nR, nP1, nP2, ls)

    # Para esta parte falta hacer los cambios de Mm/Acortamiento y comprobación anterior
    newdatafyr = regs.merge(dmF, on='Referencia', how='left')
    newdatafyr = newdatafyr.merge(dmR, on='Referencia', how='left')
    newdatafyr.columns = ['Referencia', 'FORWARD', 'REVERSE']

    newdatafyr2 = newdatafyr[(newdatafyr.FORWARD >= 1) & (newdatafyr.REVERSE >= 1)]
    lsfyr = list(newdatafyr2.Referencia)

    finallsfyr = pd.DataFrame(lsfyr, columns=['Referencia'])

    newdatafyr3 = datamerged.query('Referencia in @lsfyr')
    nF2 = newdatafyr3[newdatafyr3['Oligo'] == 'FORWARD']
    nR2 = newdatafyr3[newdatafyr3['Oligo'] == 'REVERSE']

    nF2 = nF2[['Referencia', 'Rstart', 'Acortamiento']]
    nR2 = nR2[['Referencia', 'Rstart', 'Acortamiento']]

    return newdata3, finalls, nF, nR, nP1, nP2, newdatafyr3, finallsfyr, nF2, nR2, selectedAlignments_df


def classifier(row):  # Para Inclusivity
    if (int(row.Fstart_x) < int(row.Pstart1_x) < int(row.Pstart2_x) < int(row.Rstart_y)) or \
            (int(row.Rstart_y) < int(row.Pstart2_x) < int(row.Pstart1_x) < int(row.Fstart_x)):
        row['ProbeOk'] = 'Yes'

    else:
        row['ProbeOk'] = 'No'

    return row


def primers_4_with_2_pr(filename, oligos, excel_name):
    warnings.filterwarnings('ignore')

    df = pd.read_csv(filename, sep="\t", header=None)  # Archivo .txt, resultado del blast
    df.columns = ['Oligo', 'Referencia', 'Rstart', 'Rend', 'Ostart', 'Oend', 'Mismatch', 'Kingdom', 'Group', 'Group2',
                  'Strand', 'Evalue', 'QuerySeq', 'SubjectSeq']
    df['maxL'] = df.Oend - df.Ostart + 1  # Calculamos el tamaño para comprobar si hay acortamientos
    # sentido = pd.DataFrame([df['Oligo'], df['Referencia'], df['Strand']]).T

    # Llamamos a los primers
    records = list(SeqIO.parse(oligos, "fasta"))

    olidf = pd.DataFrame({'LENGTH': [len(records[0].seq), len(records[1].seq), len(records[2].seq), len(records[3].seq)],
                          'Oligo': ['FORWARD', 'REVERSE', 'PROBE1', 'PROBE2']})

    datamerged = df.merge(olidf, on='Oligo')
    datamerged = datamerged.drop_duplicates(subset=['Oligo', 'Referencia'])
    datamerged = datamerged.reset_index()
    datamerged['Acortamiento'] = datamerged.LENGTH - datamerged.maxL

    newdata3, finalls, nF, nR, nP1, nP2, newdatafyr3, finallsfyr, nF2, nR2, selectedAlignments_df = process(datamerged)

    nP1.columns = ['Referencia', 'Pstart1', 'P-Acortamiento1', 'P-Mismatch1']
    nP2.columns = ['Referencia', 'Pstart2', 'P-Acortamiento2', 'P-Mismatch2']
    nR.columns = ['Referencia', 'Rstart', 'R-Acortamiento', 'R-Mismatch']

    selectedAlignments_df.columns = ['Referencia', 'Fstart', 'F-Acortamiento', 'F-Mismatch', 'Rstart', 'Pstart']

    # Dividimos las referencias en los diferentes primers (Forward, Reverse y Probe)
    nFF = newdata3[newdata3['Oligo'] == 'FORWARD']
    nRR = newdata3[newdata3['Oligo'] == 'REVERSE']
    nPP1 = newdata3[newdata3['Oligo'] == 'PROBE1']
    nPP2 = newdata3[newdata3['Oligo'] == 'PROBE2']

    nFF = nFF.rename(columns={'Rstart': 'Fstart'})
    nPP1 = nPP1.rename(columns={'Rstart': 'Pstart1'})
    nPP2 = nPP2.rename(columns={'Rstart': 'Pstart2'})

    # Esta tabla corresponde a la columa aling en la pestaña raw, el objetivo es comprobar la homologia de los primers
    # con el resultado obtenido despues de hacer el blast

    total = []
    tamano = []
    counter = 0

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
    for i, j, b in zip(datamerged['QuerySeq'], datamerged['Oligo'], tqdm(range(len(datamerged)))):
        replace = []
        z = 0
        aux = 0
        seq = ""

        # Comprobar homología primer a primer
        for x, y in zip(i, datamerged['SubjectSeq'][counter]):
            sentence = ""

            if len(i) == len(datamerged['SubjectSeq'][counter]):

                if x == y:  # Si el nn coincide se añade un -, representa que es homologo
                    replace.append("-")

                elif x != y:  # Si no coincide añade la letra que corresponde y suma uno para calcular el porcentaje
                    replace.append(y)
                    aux += 1

                for w in replace:  # Para que en el excel quede representado sin , o [ ]
                    sentence += str(w) + ""
                    seq = sentence

        total.append(seq)
        counter = counter + 1

        if j == 'FORWARD':  # Se calcula el porcentaje de homología para el Forward
            z = (aux / len(records[0].seq)) * 100
            z = 100 - z

        elif j == 'REVERSE':  # Se calcula el porcentaje de homología para el Reverse
            z = (aux / len(records[1].seq)) * 100
            z = 100 - z

        elif j == 'PROBE1':  # Se calcula el porcentaje de homología para el Probe
            z = (aux / len(records[2].seq)) * 100
            z = 100 - z

        elif j == 'PROBE2':  # Se calcula el porcentaje de homología para el Probe
            z = (aux / len(records[3].seq)) * 100
            z = 100 - z

        tamano.append(z)

    # Se añaden las dos columnas obtenidas (Alineamiento y porcentaje de homologia) a la tabla total
    datamerged['Aling'] = total
    datamerged['Accuracy'] = tamano

    # Evitar alineamientos multiples
    datamerged = datamerged.drop_duplicates(subset=['Oligo', 'Referencia'])

    # Esta tabla corresponde a la pestaña data
    final = selectedAlignments_df.merge(nP1[['Referencia', 'Pstart1', 'P-Acortamiento1', 'P-Mismatch1']], on='Referencia')
    final = final.merge(nP2[['Referencia', 'Pstart2', 'P-Acortamiento2', 'P-Mismatch2']], on='Referencia')
    final = final.merge(nR[['Referencia', 'Rstart', 'R-Acortamiento', 'R-Mismatch']], on='Referencia')
    final = final.drop("Rstart_x", axis=1)  # Eliminar columnas innecesarias
    final = final.drop("Pstart", axis=1)  # Eliminar columnas innecesarias
    final.rename(columns={"Rstart_x": "Rstart"}, inplace=True)

    final = final.merge(nFF[['Referencia', 'Fstart', 'Group', 'Kingdom', 'Group2']], on=['Referencia'])
    final = final.merge(nPP1[['Referencia', 'Pstart1', 'Group', 'Kingdom', 'Group2']], on=['Referencia'])
    final = final.merge(nPP2[['Referencia', 'Pstart2', 'Group', 'Kingdom', 'Group2']], on=['Referencia'])
    final = final.merge(nRR[['Referencia', 'Rstart', 'Group', 'Kingdom', 'Group2']], on=['Referencia'])
    final = final.apply(classifier, axis='columns')

    final['resta'] = final.Rstart - final.Fstart_x  # Resta para determinar si hay acortamiento
    final = final.fillna('Indeterminado')
    final = final.drop("Fstart_y", axis=1)
    final = final.drop("Pstart1_y", axis=1)
    final = final.drop("Pstart2_y", axis=1)
    final = final.drop("Group_y", axis=1)
    final = final.drop("Kingdom_y", axis=1)
    final = final.drop("Group2_y", axis=1)
    final = final.drop("Rstart", axis=1)
    final = final.loc[:, ~final.columns.duplicated()]

    # Esta tabla corresponde a la pestaña data_filter
    final_filter = final

    # Aplicar filtros para quitar aquellas secuencias que no nos interesan
    # Tamaño del amplicon de 500 pb maximo y los primers en orden correcto
    final_filter = final_filter.loc[(final_filter.ProbeOk == 'Yes') & (final_filter.resta < 500) &
                                    (final_filter.resta > -500)]

    # Eliminar duplicados para evitar alineamientos multiples
    final_filter = final_filter.drop_duplicates()

    # Corresponde a la pestaña referencias_repetidas
    referencias_repetidas = []
    localizador = []
    counter = []

    # En el caso de que haya alineamientos multiples, calcula el número de veces que se repite dicho alineamiento
    for i, a, w in zip(final_filter.Referencia, final_filter.Group2_x, tqdm(range(len(final_filter.Referencia)))):
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
    t = []

    # Separar la lista por primers
    f = list(Accuracy['Referencia'].loc[(Accuracy.Oligo == 'FORWARD')])
    r = list(Accuracy['Referencia'].loc[(Accuracy.Oligo == 'REVERSE')])
    p1 = list(Accuracy['Referencia'].loc[(Accuracy.Oligo == 'PROBE1')])
    p2 = list(Accuracy['Referencia'].loc[(Accuracy.Oligo == 'PROBE2')])
    # Tener solo en cuenta las referencias que se encuentran en las tres listas, ya que puede haber referencias a
    # secuencias que solo se encuentran en uno de los tres primers. Ese caso no amplificaria, por lo que no se tiene
    # en cuenta
    for i, b in zip(f, tqdm(range(len(f)))):
        valor = r.__contains__(i) and p1.__contains__(i) and p2.__contains__(i)

        if valor:  # Si la referencia se encuentra en los tres primers lo añade a la lista
            t.append(i)

    Acc = Accuracy
    index = []

    # Eliminar aquellas referencias que no cumplen la condicion anterior
    for i, j, b in zip(Acc['Referencia'], range(len(Acc)), tqdm(range(len(Acc['Referencia'])))):
        valor = i in t

        if not valor:
            index.append(Acc.index[j])

    Acc = Acc.drop(index)

    # Eliminar las columnas que no necesitamos en esta tabla
    Ax = Acc.drop(columns=['Rstart', 'Rend', 'Ostart', 'Oend', 'Mismatch', 'Strand', 'Evalue', 'QuerySeq', 'SubjectSeq',
                           'maxL', 'LENGTH', 'Acortamiento', 'Accuracy'])
    Ax = Ax.sort_values(by='Referencia')

    new_colum = []

    # Calcula el porcentaje para cada una de las referencias que cumplen la condicion
    for ref, c in zip(Ax['Referencia'], tqdm(range(len(Ax)))):
        suma = 0
        tam = 0

        for j, s, lens in zip(Ax['Referencia'], Ax['Aling'], Ax['Oligo']):
            if ref == j:
                aux = s.count('-')  # Cuenta el numero de nn homologos
                suma = suma + aux  # Como son 3 primers tiene que tener en cuenta los tres

                if lens == "FORWARD":  # Suma el tamaño del primer Forward
                    tam = tam + len(records[0].seq)

                if lens == "REVERSE":  # Suma el tamaño del primer Reverse
                    tam = tam + len(records[1].seq)

                if lens == "PROBE1":  # Suma el tamaño del primer Probe
                    tam = tam + len(records[2].seq)

                if lens == "PROBE2":  # Suma el tamaño del primer Probe
                    tam = tam + len(records[3].seq)

        total = (suma / tam) * 100  # Calcula el porcentaje
        new_colum.append(total)

    Ax['Total_Accuracy'] = new_colum  # Añade el porcentaje a la tabla

    Ax_final = Ax.drop_duplicates(subset=['Referencia'])
    Ax_final = Ax_final.drop(columns=['Oligo', 'Aling'])

    # Esta tabla corresponde a la pestaña resume
    final2 = final.drop_duplicates(subset=['Referencia'])
    final2 = final2.groupby(['Kingdom_x', 'Group_x', 'Group2_x']).Referencia.count()

    # Esta tabla corresponde a la pestaña ProbeOk
    final3 = final.drop_duplicates(subset=['Referencia'])
    final3 = final3[(final3.ProbeOk == 'Yes')].groupby(['Kingdom_x', 'Group_x', 'Group2_x']).Referencia.count()

    # Esta tabla corresponde a la pestaña mm
    final3mm_filter = final
    # Aplicamos los filtros del tamaño del amplicon y tener los primers en el orden correcto
    final3mm_filter = final3mm_filter.loc[(final3mm_filter.ProbeOk == 'Yes') &
                                          (final3mm_filter.resta < 500) &
                                          (final3mm_filter.resta > -500)]
    # Contamos el numero de veces que se repite cada alineamiento
    final3mm = final3mm_filter.groupby(['Kingdom_x', 'Group_x', 'Group2_x', 'F-Acortamiento', 'F-Mismatch',
                                        'P-Acortamiento1', 'P-Mismatch1', 'P-Acortamiento2', 'P-Mismatch2',
                                        'R-Acortamiento', 'R-Mismatch', 'resta']).Referencia.count()

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
