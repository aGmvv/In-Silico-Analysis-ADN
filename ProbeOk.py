"""
coding=utf-8
author = Alejandro González

Estas funciones comprueban si la sonda se encuentra entre los primers forward y reverse.
Esto lo revisa mirando las posiciones de inicio y finas, que la sonda se encuentre entre los primers forward y
reverse es indispensable para la amplificación de la PCR.
"""


def classifier(row):  # Para Inclusivity
    if (int(row.Fstart_x) < int(row.Pstart_y) < int(row.Rstart_y)) or (int(row.Rstart_y)
                                                                       < int(row.Pstart_y) < int(row.Fstart_x)):
        row['ProbeOk'] = 'Yes'

    else:
        row['ProbeOk'] = 'No'

    return row


def classifier2(row):  # Para Cross Reactivity
    if (int(row.Fstart) < int(row.Pstart_y) < int(row.Rstart_y)) or (int(row.Rstart_y)
                                                                     < int(row.Pstart_y) < int(row.Fstart)):
        row['ProbeOk'] = 'Yes'

    else:
        row['ProbeOk'] = 'No'

    return row
