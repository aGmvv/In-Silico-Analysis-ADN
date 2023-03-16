import pandas as pd


def ID_primers():
    data = pd.read_excel('Z:/MDx/qPCR Proveedores/OLIGOS/Lista secuencias definitivas.xlsx',
                         sheet_name='Primers patógenos', skiprows=1)
    data['Ref-COD 3 LETRAS'].fillna(method='ffill', inplace=True)

    with open('Z:/MDx/qPCR Proveedores/OLIGOS/Lista primers.txt', 'w') as archivo:
        for ref, f, r, p in zip(data['Ref-COD 3 LETRAS'], data['Primer Fw (5´-3´)'],
                                data['Reverse (5´-3´)'], data['Sonda ']):
            if str(f) != 'nan':
                archivo.write('>' + str(ref) + '-f' + '\n')
                archivo.write(str(f) + '\n' + '\n')

            if str(r) != 'nan':
                archivo.write('>' + str(ref) + '-r' + '\n')
                archivo.write(str(r) + '\n' + '\n')

            if str(p) != 'nan':
                archivo.write('>' + str(ref) + '-p' + '\n')
                archivo.write(str(p) + '\n' + '\n')

            archivo.write('----------------------------------------------------------------------------------'
                          + '\n' + '\n')
