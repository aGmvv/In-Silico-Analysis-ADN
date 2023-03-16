import time
import pandas as pd


def bucle(colum):

    nuevo_df = []
    for col in (colum['Query']):
        nuevo_df.append(list(col))

    df = pd.DataFrame(nuevo_df)

    resultados = []

    for col2 in df.columns:
        # Contar la frecuencia de cada valor único en la columna
        frecuencia = df[col2].value_counts()
        # Obtener el valor más frecuente y su frecuencia
        letra_mas_frecuente = frecuencia.index[0]
        frecuencia_letra_mas_frecuente = frecuencia.iloc[0]
        # Agregar el resultado a la lista de resultados
        resultados.append((letra_mas_frecuente, frecuencia_letra_mas_frecuente))

    new_df = pd.DataFrame(resultados)
    new_df.columns = ['Nucleotido', 'Cantidad de veces']
    new_df = new_df.T

    return new_df


def contar(file):

    tic = time.perf_counter()

    g1 = pd.read_excel(file, sheet_name='Bad')
    f = g1[g1['Oliname'] == 'Forward']
    f = f[f['MultiClass'] != '4']
    r = g1[g1['Oliname'] == 'Reverse']
    r = r[r['MultiClass'] != '4']
    p = g1[g1['Oliname'] == 'Probe']
    p = p[p['MultiClass'] != '4']

    f = f.loc[:, ['Query']]
    r = r.loc[:, ['Query']]
    p = p.loc[:, ['Query']]

    f = f.dropna(subset=['Query'])
    r = r.dropna(subset=['Query'])
    p = p.dropna(subset=['Query'])

    f_excel = bucle(f)
    r_excel = bucle(r)
    p_excel = bucle(p)

    f_excel.to_excel('archivo_creado.xlsx', sheet_name='forward')
    with pd.ExcelWriter('archivo_creado.xlsx', mode='a', engine="openpyxl") as writer:
        r_excel.to_excel(writer, sheet_name='reverse')
        p_excel.to_excel(writer, sheet_name='probe')

    print("FINISH")
    toc = time.perf_counter()
    print(f"Code in {toc - tic:0.4f} seconds")
