"""
coding=utf-8
author = Alejandro González
"""

import pandas as pd
from aligmentselection import aligmentselection, aligmentselection2


def process(datamerged):
    # Seleccionar las filas de datos que corresponden a Oligo "FORWARD" y agrupar por Referencia
    interF = datamerged.query('Oligo == "FORWARD"')
    dmF = interF.groupby(['Referencia']).Oligo.count()
    dmF = pd.DataFrame(dmF)  # Convertir el resultado en un DataFrame
    dmF = dmF.reset_index()  # Reiniciar el índice del DataFrame

    # Seleccionar las filas de datos que corresponden a Oligo "REVERSE" y agrupar por Referencia
    interR = datamerged.query('Oligo == "REVERSE"')
    dmR = interR.groupby(['Referencia']).Oligo.count()
    dmR = pd.DataFrame(dmR)  # Convertir el resultado en un DataFrame
    dmR = dmR.reset_index()  # Reiniciar el índice del DataFrame

    # Seleccionar las filas de datos que corresponden a Oligo "PROBE" y agrupar por Referencia
    interP = datamerged.query('Oligo == "PROBE"')
    dmP = interP.groupby(['Referencia']).Oligo.count()
    dmP = pd.DataFrame(dmP)  # Convertir el resultado en un DataFrame
    dmP = dmP.reset_index()  # Reiniciar el índice del DataFrame

    # Obtener las referencias únicas y crear un nuevo DataFrame con ellas
    regs = datamerged.Referencia.unique()
    regs = pd.DataFrame(regs, columns=['Referencia'])

    # Unir los DataFrames dmF, dmR, y dmP en un nuevo DataFrame "newdata"
    newdata = regs.merge(dmF, on='Referencia', how='left')
    newdata = newdata.merge(dmR, on='Referencia', how='left')
    newdata = newdata.merge(dmP, on='Referencia', how='left')
    newdata.columns = ['Referencia', 'FORWARD', 'REVERSE', 'PROBE']

    # Seleccionar las referencias que tengan al menos una entrada para cada tipo de oligo
    newdata2 = newdata[(newdata.FORWARD >= 1) & (newdata.REVERSE >= 1) & (newdata.PROBE >= 1)]
    ls = list(newdata2.Referencia)  # Obtener una lista de referencias seleccionadas

    # Se crea un dataframe 'finalls' con la columna 'Referencia' y los valores de la lista 'ls'
    finalls = pd.DataFrame(ls, columns=['Referencia'])
    # Se filtra el dataframe 'datamerged' por las referencias en 'ls' y se guarda en 'newdata3'
    newdata3 = datamerged.query('Referencia in @ls')

    # Se separan las filas del dataframe 'newdata3' correspondientes a las columnas 'Oligo' 'FORWARD', 'REVERSE',
    # y 'PROBE'
    nF = newdata3[newdata3['Oligo'] == 'FORWARD']
    nR = newdata3[newdata3['Oligo'] == 'REVERSE']
    nP = newdata3[newdata3['Oligo'] == 'PROBE']

    # Se seleccionan las columnas 'Referencia', 'Rstart', 'Acortamiento', 'Mismatch' de cada uno de los
    # dataframes separados
    nF = nF[['Referencia', 'Rstart', 'Acortamiento', 'Mismatch']]
    nR = nR[['Referencia', 'Rstart', 'Acortamiento', 'Mismatch']]
    nP = nP[['Referencia', 'Rstart', 'Acortamiento', 'Mismatch']]

    #  Se llama a la función 'aligmentselection' y se guarda el resultado en 'selectedAlignments_df'
    selectedAlignments_df = aligmentselection(nF, nR, nP, ls)

    # Se hace un merge entre los dataframes 'regs' y 'dmF' por la columna 'Referencia' y se guarda en 'newdatafyr'
    newdatafyr = regs.merge(dmF, on='Referencia', how='left')
    # Se hace un merge entre el dataframe 'newdatafyr' y 'dmR' por la columna 'Referencia' y se guarda en 'newdatafyr'
    newdatafyr = newdatafyr.merge(dmR, on='Referencia', how='left')
    # Se renombran las columnas 'FORWARD' y 'REVERSE'
    newdatafyr.columns = ['Referencia', 'FORWARD', 'REVERSE']

    # Se filtra el dataframe 'newdatafyr' por las referencias que tienen valores mayores o iguales a 1 en las
    # columnas 'FORWARD' y 'REVERSE'
    newdatafyr2 = newdatafyr[(newdatafyr.FORWARD >= 1) & (newdatafyr.REVERSE >= 1)]
    # Se obtiene una lista 'lsfyr' con las referencias filtradas
    lsfyr = list(newdatafyr2.Referencia)

    # Se crea un dataframe 'finallsfyr' con la columna 'Referencia' y los valores de la lista 'lsfyr'
    finallsfyr = pd.DataFrame(lsfyr, columns=['Referencia'])

    # Se filtran los datos en "datamerged" para aquellos cuyo valor de "Referencia" se encuentre en "lsfyr"
    # y se almacenan en "newdatafyr3".
    newdatafyr3 = datamerged.query('Referencia in @lsfyr')
    # Se separan los datos de la columna "Oligo" para aquellos que tengan el valor "FORWARD" y se almacenan en "nF2".
    nF2 = newdatafyr3[newdatafyr3['Oligo'] == 'FORWARD']
    # Se separan los datos de la columna "Oligo" para aquellos que tengan el valor "REVERSE" y se almacenan en "nR2".
    nR2 = newdatafyr3[newdatafyr3['Oligo'] == 'REVERSE']

    # Se seleccionan solo las columnas relevantes para "nF2" y "nR2".
    nF2 = nF2[['Referencia', 'Rstart', 'Acortamiento']]
    nR2 = nR2[['Referencia', 'Rstart', 'Acortamiento']]

    return newdata3, finalls, nF, nR, nP, newdatafyr3, finallsfyr, nF2, nR2, selectedAlignments_df


def process2(datamerged):
    # Seleccionar las filas de datos que corresponden a Oligo "FORWARD" y agrupar por Referencia
    interF = datamerged.query('Oliname == "Forward"')
    dmF = interF.groupby(['Registro']).Oliname.count()
    dmF = pd.DataFrame(dmF)  # Convertir el resultado en un DataFrame
    dmF = dmF.reset_index()  # Reiniciar el índice del DataFrame

    # Seleccionar las filas de datos que corresponden a Oligo "REVERSE" y agrupar por Referencia
    interR = datamerged.query('Oliname == "Reverse"')
    dmR = interR.groupby(['Registro']).Oliname.count()
    dmR = pd.DataFrame(dmR)  # Convertir el resultado en un DataFrame
    dmR = dmR.reset_index()  # Reiniciar el índice del DataFrame

    # Seleccionar las filas de datos que corresponden a Oligo "PROBE" y agrupar por Referencia
    interP = datamerged.query('Oliname == "Probe"')
    dmP = interP.groupby(['Registro']).Oliname.count()
    dmP = pd.DataFrame(dmP)  # Convertir el resultado en un DataFrame
    dmP = dmP.reset_index()  # Reiniciar el índice del DataFrame

    # Obtener las referencias únicas y crear un nuevo DataFrame con ellas
    regs = datamerged.Registro.unique()
    regs = pd.DataFrame(regs, columns=['Registro'])

    # Unir los DataFrames dmF, dmR, y dmP en un nuevo DataFrame "newdata"
    newdata = regs.merge(dmF, on='Registro', how='left')
    newdata = newdata.merge(dmR, on='Registro', how='left')
    newdata = newdata.merge(dmP, on='Registro', how='left')
    newdata.columns = ['Registro', 'Forward', 'Reverse', 'Probe']

    # Seleccionar las referencias que tengan al menos una entrada para cada tipo de oligo
    newdata2 = newdata[(newdata.Forward >= 1) & (newdata.Reverse >= 1) & (newdata.Probe >= 1)]
    ls = list(newdata2.Registro)  # Obtener una lista de referencias seleccionadas

    # Se crea un dataframe 'finalls' con la columna 'Referencia' y los valores de la lista 'ls'
    finalls = pd.DataFrame(ls, columns=['Registro'])
    # Se filtra el dataframe 'datamerged' por las referencias en 'ls' y se guarda en 'newdata3'
    newdata3 = datamerged.query('Registro in @ls')

    # Se separan las filas del dataframe 'newdata3' correspondientes a las columnas 'Oligo' 'Forward', 'Reverse',
    # y 'Probe'
    nF = newdata3[newdata3['Oliname'] == 'Forward']
    nR = newdata3[newdata3['Oliname'] == 'Reverse']
    nP = newdata3[newdata3['Oliname'] == 'Probe']

    # Se seleccionan las columnas 'Referencia', 'Rstart', 'Acortamiento', 'Mismatch' de cada uno de los
    # dataframes separados
    nF = nF[['Registro', 'Inicio', 'Acortamiento', 'MM']]
    nR = nR[['Registro', 'Inicio', 'Acortamiento', 'MM']]
    nP = nP[['Registro', 'Inicio', 'Acortamiento', 'MM']]

    #  Se llama a la función 'aligmentselection' y se guarda el resultado en 'selectedAlignments_df'
    selectedAlignments_df = aligmentselection2(nF, nR, nP, ls)

    # Se hace un merge entre los dataframes 'regs' y 'dmF' por la columna 'Referencia' y se guarda en 'newdatafyr'
    newdatafyr = regs.merge(dmF, on='Registro', how='left')
    # Se hace un merge entre el dataframe 'newdatafyr' y 'dmR' por la columna 'Referencia' y se guarda en 'newdatafyr'
    newdatafyr = newdatafyr.merge(dmR, on='Registro', how='left')
    # Se renombran las columnas 'Forward' y 'Reverse'
    newdatafyr.columns = ['Registro', 'Forward', 'Reverse']

    # Se filtra el dataframe 'newdatafyr' por las referencias que tienen valores mayores o iguales a 1 en las
    # columnas 'Forward' y 'Reverse'
    newdatafyr2 = newdatafyr[(newdatafyr.Forward >= 1) & (newdatafyr.Reverse >= 1)]
    # Se obtiene una lista 'lsfyr' con las referencias filtradas
    lsfyr = list(newdatafyr2.Registro)

    # Se crea un dataframe 'finallsfyr' con la columna 'Referencia' y los valores de la lista 'lsfyr'
    finallsfyr = pd.DataFrame(lsfyr, columns=['Registro'])

    # Se filtran los datos en "datamerged" para aquellos cuyo valor de "Referencia" se encuentre en "lsfyr"
    # y se almacenan en "newdatafyr3".
    newdatafyr3 = datamerged.query('Registro in @lsfyr')
    # Se separan los datos de la columna "Oligo" para aquellos que tengan el valor "FORWARD" y se almacenan en "nF2".
    nF2 = newdatafyr3[newdatafyr3['Oliname'] == 'Forward']
    # Se separan los datos de la columna "Oligo" para aquellos que tengan el valor "REVERSE" y se almacenan en "nR2".
    nR2 = newdatafyr3[newdatafyr3['Oliname'] == 'Reverse']

    # Se seleccionan solo las columnas relevantes para "nF2" y "nR2".
    nF2 = nF2[['Registro', 'Inicio', 'Acortamiento']]
    nR2 = nR2[['Registro', 'Inicio', 'Acortamiento']]

    return newdata3, finalls, nF, nR, nP, newdatafyr3, finallsfyr, nF2, nR2, selectedAlignments_df
