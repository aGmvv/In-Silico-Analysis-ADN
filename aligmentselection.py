"""
coding=utf-8
author = Alejandro González
"""

import numpy as np
import pandas as pd
from tqdm import tqdm


def aligmentselection(Fw_df, Rv_df, Pr_df, ls):
    # Crear un nuevo dataframe a partir de Fw_df con los nombres de columna renombrados
    new_df = pd.DataFrame(Fw_df)
    new_df = new_df.rename(columns={'Rstart': 'Fstart'})

    #  A partir de la lista con las referencias finales, se crean diccionarios de la siguiente manera:
    dR = {}  # dR['Referencia'] = (lista de posibles alineamientos → todas las posiciones de inicio que
    #  tiene el Rv en cada secuencia)
    for ref, c in zip(ls, tqdm(range(len(ls)))):
        dR[ref] = list(Rv_df['Rstart'].loc[Rv_df['Referencia'] == ref])

    dP = {}  # dP['Referencia'] = (lista de posibles alineamientos → todas las posiciones de inicio que
    #  tiene la Pr en cada secuencia)
    for ref, c in zip(ls, tqdm(range(len(ls)))):
        dP[ref] = list(Pr_df['Rstart'].loc[Pr_df['Referencia'] == ref])
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
            PrStart_temp = [Rstart for Rstart in dP[ref] if (fw < Rstart < rv) or (rv < Rstart < fw)]

            if len(PrStart_temp) >= 1:
                new_df['Pstart'].loc[new_df['Fstart'] == fw] = PrStart_temp[0]

    new_df = new_df.dropna()

    return new_df


def aligmentselection2(Fw_df, Rv_df, Pr_df, ls):
    """
    Esta función tiene como objetivo seleccionar los alineamientos apropiados entre tres dataframes: Fw_df, Rv_df y
    Pr_df, donde cada dataframe contiene información de una dirección de lectura de una secuencia de ADN. La lista
    "ls" es una lista de referencias que se utilizarán para seleccionar los alineamientos apropiados.

    La función comienza creando un nuevo dataframe a partir de Fw_df y renombrando la columna "Inicio" como "Fstart".
    Luego, crea dos diccionarios, dR y dP, que contienen información sobre las posibles posiciones de inicio de las
    lecturas Rv y Pr en cada secuencia.

    A continuación, se realiza una comprobación de que el amplicón entre la lectura Fw y Rv tenga el tamaño adecuado.
    Se recorre la lista de referencias y, para cada una, se recorre la lista de posibles alineamientos de la lectura Fw.
    Si existe al menos un alineamiento válido entre la lectura Fw y Rv, se elige el primero y se añade su posición
    inicial a una nueva columna en el dataframe. Finalmente, se eliminan las filas con valores NaN en el dataframe.

    Por último, se realiza una comprobación de que la sonda esté en medio entre las lecturas Fw y Rv. Se recorre
    nuevamente la lista de referencias y, para cada una, se crea una tupla con las posiciones iniciales de las lecturas
    Fw y Rv. Luego, se recorre esta tupla y se compara con la lista de posibles alineamientos de la lectura Pr.
    Si existe al menos un alineamiento válido entre la sonda y Fw/Rv, se elige el primero y se añade su posición
    inicial a una nueva columna en el dataframe. Finalmente, se eliminan las filas con valores NaN en el dataframe y
    se devuelve el dataframe resultante.
    :param Fw_df:
    :param Rv_df:
    :param Pr_df:
    :param ls:
    :return:
    """
    # Crea un nuevo dataframe a partir de Fw_df
    new_df = pd.DataFrame(Fw_df)
    # Renombra la columna 'Inicio' a 'Fstart'
    new_df = new_df.rename(columns={'Inicio': 'Fstart'})
    #  A partir de la lista con las referencias finales, se crean diccionarios de la siguiente manera:
    dR = {}  # dR['Referencia'] = (lista de posibles alineamientos → todas las posiciones de inicio que
    #  tiene el Rv en cada secuencia)

    for ref, c in zip(ls, tqdm(range(len(ls)))):
        dR[ref] = list(Rv_df['Inicio'].loc[Rv_df['Registro'] == ref])
    dP = {}  # dP['Referencia'] = (lista de posibles alineamientos → todas las posiciones de inicio que
    #  tiene la Pr en cada secuencia)

    for ref, c in zip(ls, tqdm(range(len(ls)))):
        dP[ref] = list(Pr_df['Inicio'].loc[Pr_df['Registro'] == ref])
    #  Comprobación de que el amplicon entre Fw y Rv tiene el tamaño adecuado
    new_df['Inicio'] = np.nan

    for ref, c in zip(ls, tqdm(range(len(ls)))):  # Se recorre la lista de reference
        Fstarts = list(new_df['Fstart'].loc[new_df['Registro'] == ref])

        for st in Fstarts:  # Para cada referencia, se recorre la si lista de posibles alineamientos del forward
            RvStart_temp = [Rstart for Rstart in dR[ref] if abs(Rstart - st) < 500]

            if len(RvStart_temp) >= 1:  # Si para un Fw hay varios Rv posibles, se elige el primero
                new_df['Inicio'].loc[new_df['Fstart'] == st] = RvStart_temp[0]  # Se añade la posicion start del
                #  Fw a una nueva col

    new_df = new_df.dropna()
    #  Comprobación de que la sonda está en medio entre Fw y Rv
    new_df['Pstart'] = np.nan

    for ref, c in zip(ls, tqdm(range(len(ls)))):
        Fstarts = list(new_df['Fstart'].loc[new_df['Registro'] == ref])
        RVstarts = list(new_df['Inicio'].loc[new_df['Registro'] == ref])
        tuple_FwRv = tuple(zip(Fstarts, RVstarts))

        for tupl in tuple_FwRv:
            fw = tupl[0]
            rv = tupl[1]
            PrStart_temp = [Rstart for Rstart in dP[ref] if (fw < Rstart < rv) or (rv < Rstart < fw)]

            if len(PrStart_temp) >= 1:
                new_df['Pstart'].loc[new_df['Fstart'] == fw] = PrStart_temp[0]

    new_df = new_df.dropna()
    return new_df
