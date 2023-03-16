import warnings
import pandas as pd
from tqdm import tqdm


def mutation_differentation(filename, excel_name):
    warnings.filterwarnings('ignore')
    file = pd.read_excel(filename, sheet_name='Bad')

    data = file.copy()
    data = data[data["MultiClass"] != 4]
    data = data.drop(['Ensayo', 'Score', 'MM', 'Inicio', 'Fin', 'Query', 'Variant', 'Chain', '¿Multiple?',
                      'CHAIN2', 'SeqGen', 'MultiClass',
                      'ComboVariant', 'SingleClass', 'Max E1',
                      'Max -', 'Max -.1', 'Max -.2'], axis=1)
    registro = []
    fw = []
    rv = []
    pr = []
    total = len(data) / 3

    for i, c in zip(range(0, len(data), 3), tqdm(range(int(total)))):
        group = data[i:i + 3]
        group.reset_index(inplace=True)
        registro.append(str(group['Registro'][0]))
        fw.append(str(group['Result'][0]))
        rv.append(str(group['Result'][1]))
        pr.append(str(group['Result'][2]))

    final = pd.DataFrame(columns=['Registro', 'Fw', 'Rv', 'Pr'])
    final['Registro'] = registro
    final['Fw'] = fw
    final['Rv'] = rv
    final['Pr'] = pr

    result = final.groupby(['Fw', 'Rv', 'Pr']).size().reset_index(name='counts')

    print("Creando el excel ...")
    # Crear Excel
    final.to_excel(excel_name, sheet_name='all')
    with pd.ExcelWriter(excel_name, mode='a', engine="openpyxl") as writer:
        result.to_excel(writer, sheet_name='counts')

    print("Ha acabado")
