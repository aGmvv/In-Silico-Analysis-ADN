"""
coding=utf-8
author = Alejandro González
"""

import os
import time
import warnings
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from filtros import process2
from ProbeOk import classifier2


def porcentaje_homologia(registro, result, oliname, newcolum, records, final_filter1):
    x = pd.DataFrame(registro).drop_duplicates()
    for ref, c in zip(x['Registro'], tqdm(range(len(x)))):  # Calcula el porcentaje de homologia

        if ref in list(final_filter1.Registro):
            suma = 0
            tam = 0
            for j, s, lens in zip(registro, result, oliname):  # Compara cada letra de los primers

                if ref == j:
                    aux = str(s).count('-')  # Si es la misma es homologo
                    suma = suma + int(aux)

                    if lens == "Forward":  # Calcula la distancia para Forward
                        tam = tam + len(records[0].seq)

                    if lens == "Reverse":  # Calcula la distancia para Reverse
                        tam = tam + len(records[1].seq)

                    if lens == "Probe":  # Calcula la distancia para Probe
                        tam = tam + len(records[2].seq)

            total = (suma / tam) * 100  # Calcula el porcentaje total
            for w in range(3):
                newcolum.append(total)
        else:
            for w in range(3):
                newcolum.append('No ok')
    return newcolum


def preparar_cross(file_name, excel_name1, excel_name2, excel_name3, excel_name4,
                   oligos1, oligos2, oligos3, oligos4, diana):
    warnings.filterwarnings('ignore')

    read_file = pd.read_excel(file_name, sheet_name='Bad')  # Lee el excel y crea una copia en txt para manejarlo
    read_file.to_csv(r'file.txt', index=False, sep="\t", encoding="utf-8")
    df = pd.read_csv(r'file.txt', sep="\t", encoding="utf-8")
    os.remove('file.txt')

    if diana == 1:  # Cuando solo hay una diana
        e1 = df.loc[(df.Ensayo == 'E1')]

        # Elimina los duplicados y todas las columnas que no son necesarias
        e1 = e1.drop_duplicates(subset=['Registro', 'Oliname'])

        e1 = e1.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max -', 'Max -.1',
                              'Max -.2'])

        final_cross(excel_name1, oligos1, e1)

    if diana == 2:  # Para dos dianas
        e1 = df.loc[(df.Ensayo == 'E1')]
        e2 = df.loc[(df.Ensayo == 'E2')]

        # Elimina los duplicados y todas las columnas que no son necesarias
        e1 = e1.drop_duplicates(subset=['Registro', 'Oliname'])
        e2 = e2.drop_duplicates(subset=['Registro', 'Oliname'])

        e1 = e1.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max E2', 'Max -',
                              'Max -.1'])
        e2 = e2.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max E2', 'Max -',
                              'Max -.1'])

        # Hacer el codigo en paralelo
        final_cross(excel_name1, oligos1, e1)
        final_cross(excel_name2, oligos2, e2)

    if diana == 3:  # Para tres dianas
        e1 = df.loc[(df.Ensayo == 'E1')]
        e2 = df.loc[(df.Ensayo == 'E2')]
        e3 = df.loc[(df.Ensayo == 'E3')]

        # Elimina los duplicados y todas las columnas que no son necesarias
        e1 = e1.drop_duplicates(subset=['Registro', 'Oliname'])
        e2 = e2.drop_duplicates(subset=['Registro', 'Oliname'])
        e3 = e3.drop_duplicates(subset=['Registro', 'Oliname'])

        e1 = e1.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max -', 'Max E2', 'Max E3'])
        e2 = e2.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max -', 'Max E2', 'Max E3'])
        e3 = e3.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max -', 'Max E2', 'Max E3'])

        # Hacer el codigo en paralelo
        final_cross(excel_name1, oligos1, e1)
        final_cross(excel_name2, oligos2, e2)
        final_cross(excel_name3, oligos3, e3)

    if diana == 4:
        e1 = df.loc[(df.Ensayo == 'E1')]
        e2 = df.loc[(df.Ensayo == 'E2')]
        e3 = df.loc[(df.Ensayo == 'E3')]
        e4 = df.loc[(df.Ensayo == 'E4')]

        # Elimina los duplicados y todas las columnas que no son necesarias
        e4 = e4.drop_duplicates(subset=['Registro', 'Oliname'])
        e1 = e1.drop_duplicates(subset=['Registro', 'Oliname'])
        e2 = e2.drop_duplicates(subset=['Registro', 'Oliname'])
        e3 = e3.drop_duplicates(subset=['Registro', 'Oliname'])

        e1 = e1.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max E2', 'Max E3',
                              'Max E4'])
        e2 = e2.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max E2', 'Max E3',
                              'Max E4'])
        e3 = e3.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max E2', 'Max E3',
                              'Max E4'])
        e4 = e4.drop(columns=['Ensayo', 'Score', 'Target', 'Query', '¿Multiple?', 'CHAIN2', 'SeqGen', 'MultiClass',
                              'ComboVariant', 'SingleClass', 'Variant', 'Chain', 'Max E1', 'Max E2', 'Max E3',
                              'Max E4'])

        final_cross(excel_name1, oligos1, e1)
        final_cross(excel_name2, oligos2, e2)
        final_cross(excel_name3, oligos3, e3)
        final_cross(excel_name4, oligos4, e4)


def final_cross(excel_name1, oligos1, ensayo):
    new_colum1 = []
    warnings.filterwarnings('ignore')
    tic = time.perf_counter()

    """
    Esta funcion sigue los mismos pasos que la funcion Inclusivity, pero esta adaptado al tipo de tabla que se le 
    da como archivo.
    """

    ensayo.columns = ['Registro', 'Oliname', 'MM', 'Inicio', 'Fin', 'Result']
    ensayo['maxL'] = ensayo.Fin - ensayo.Inicio + 1  # Calculamos el tamaño para comprobar si hay acortamientos

    # Llamamos a los primers
    records1 = list(SeqIO.parse(oligos1, "fasta"))
    olidf1 = pd.DataFrame({'LENGTH': [len(records1[0].seq), len(records1[1].seq), len(records1[2].seq)],
                           'Oliname': ['Forward', 'Reverse', 'Probe']})

    datamerged1 = ensayo.merge(olidf1, on='Oliname')  # Unimos los datos que hemos obtenido del blast con los primers

    #  Esta col cuenta los mismatch del principio o final mirando si coinciden longitudes
    datamerged1['Acortamiento'] = datamerged1.LENGTH - datamerged1.maxL
    for i in range(len(datamerged1['Acortamiento'])):
        if datamerged1['Acortamiento'][i] < 0:
            datamerged1['Acortamiento'][i] = 0

    newdata31, finalls1, nF1, nR1, nP1, newdatafyr31, finallsfyr1, nF21, nR21, selectedAlignments_df1 = process2(
        datamerged1)

    nP1.columns = ['Registro', 'Pstart', 'P-Acortamiento', 'P-Mismatch']
    nR1.columns = ['Registro', 'Rstart', 'R-Acortamiento', 'R-Mismatch']

    selectedAlignments_df1.columns = ['Registro', 'Fstart', 'F-Acortamiento', 'F-Mismatch', 'Rstart', 'Pstart']

    # Dividimos las referencias en los diferentes primers (Forward, Reverse y Probe)
    nFF1 = newdata31[newdata31['Oliname'] == 'Forward']
    nPP1 = newdata31[newdata31['Oliname'] == 'Probe']
    nRR1 = newdata31[newdata31['Oliname'] == 'Reverse']

    nFF1 = nFF1.rename(columns={'Inico': 'Fstart'})
    nPP1 = nPP1.rename(columns={'Inico': 'Pstart'})

    # Esta tabla corresponde a la pestaña data
    final1 = selectedAlignments_df1.merge(nP1[['Registro', 'Pstart', 'P-Acortamiento', 'P-Mismatch']], on='Registro')

    final1 = final1.merge(nR1[['Registro', 'Rstart', 'R-Acortamiento', 'R-Mismatch']], on='Registro')
    final1 = final1.drop("Rstart_x", axis=1)  # Eliminar columnas innecesarias
    final1 = final1.drop("Pstart_x", axis=1)  # Eliminar columnas innecesarias
    final1.rename(columns={"Rstart_x": "Rstart", "Pstart_x": "Pstart"}, inplace=True)

    final1 = final1.merge(nFF1[['Registro', 'Inicio']], on=['Registro'])
    final1 = final1.merge(nPP1[['Registro', 'Inicio']], on=['Registro'])
    final1 = final1.merge(nRR1[['Registro', 'Inicio']], on=['Registro'])

    final1 = final1.apply(classifier2, axis='columns')

    """
    En el caso de que ninguna secuencia cumpla los requisitos anteriores, se añade una columna vacia para que el 
    codigo pueda continuar sin dar errores
    """

    if len(final1) == 0:
        final1['ProbeOk'] = ""

    final1['resta'] = final1.Rstart_y - final1.Fstart  # Resta para determinar si hay acortamiento

    final1 = final1.fillna('Indeterminado')
    final1 = final1.drop_duplicates()

    # Esta tabla corresponde a la pestaña data_filter
    final_filter1 = final1

    # Aplicar filtros para quitar aquellas secuencias que no nos interesan
    # Tamaño del amplicon de 500 pb maximo y los primers en orden correcto
    final_filter1 = final_filter1.loc[(final_filter1.ProbeOk == 'Yes') & (final_filter1.resta < 500)
                                      & (final_filter1.resta > -500)]

    final_filter1 = final_filter1.drop_duplicates()

    Accuracy1 = datamerged1

    # Separar la lista por primers
    f1 = list(Accuracy1['Registro'].loc[(Accuracy1.Oliname == 'Forward')])
    r1 = list(Accuracy1['Registro'].loc[(Accuracy1.Oliname == 'Reverse')])
    p1 = list(Accuracy1['Registro'].loc[(Accuracy1.Oliname == 'Probe')])

    # Tener solo en cuenta las referencias que se encuentran en las tres listas, ya que puede haber referencias a
    # secuencias que solo se encuentran en uno de los tres primers. Ese caso no amplificaria, por lo que no se tiene
    # en cuenta
    t1 = [i for i in f1 if i in r1 and i in p1]

    Acc1 = Accuracy1

    # Eliminar aquellas referencias que no cumplen la condicion anterior
    index = [Acc1.index[j] for i, j in zip(Acc1['Registro'], range(len(Acc1))) if i not in t1]
    Acc1 = Acc1.drop(index)

    # Eliminar las columnas que no necesitamos en esta tabla
    Ax1 = Acc1.drop(columns=['MM', 'Inicio', 'Fin', 'maxL', 'LENGTH', 'Acortamiento'])
    Ax1 = Ax1.sort_values(by='Registro')

    # Porcentaje de Homología total
    new_colum1 = porcentaje_homologia(Ax1['Registro'], Ax1['Result'], Ax1['Oliname'], new_colum1, records1,
                                      final_filter1)
    try:
        Ax1['Total_Accuracy'] = new_colum1
        Ax_final1 = Ax1.drop_duplicates(subset=['Registro'])
        Ax_final1 = Ax_final1.drop(columns=['Oliname', 'Result'])
        datamerged1.to_excel(excel_name1, sheet_name='raw')

        with pd.ExcelWriter(excel_name1, mode='a', engine="openpyxl") as writer:
            final1.to_excel(writer, sheet_name='datos')
            final_filter1.to_excel(writer, sheet_name='datos_filter')
            Ax_final1.to_excel(writer, sheet_name='Porcentaje')
    except:

        # Crear excel con los resultados
        datamerged1.to_excel(excel_name1, sheet_name='raw')
        with pd.ExcelWriter(excel_name1, mode='a', engine="openpyxl") as writer:
            final1.to_excel(writer, sheet_name='datos')
            final_filter1.to_excel(writer, sheet_name='datos_filter')

    toc = time.perf_counter()
    print(f"Code in {toc - tic:0.4f} seconds")
    print("FINISH Cross")
