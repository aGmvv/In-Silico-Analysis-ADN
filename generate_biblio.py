import os
import itertools
import collections
from tqdm import tqdm
from PyPDF2 import PdfReader
import matplotlib.pyplot as plt
import codecs
import pandas as pd
# find "Z:\MDx\qPCR Proyectos" -type d -name "Bibliografía" -execdir sh -c 'cp -t Z:/Temporal/Alejandro/Folder/ "{}"/*.pdf' \;
# Abrir el archivo PDF
lista = []
entrada = codecs.open("export.xml", "r", encoding="utf-8")
# for i in entrada:
#
#     linia = i.rstrip().lstrip()
#     if linia.startswith('<periodical><full-title>'):
#         x = linia.replace('<periodical><full-title>', '')
#         x = x.replace('</full-title></periodical>', '')
#         lista.append(x)
#
# counter = collections.Counter(lista)
# dicc = dict(sorted(counter.items(), key=lambda item: item[1], reverse=True))
# # diccionario = dict(itertools.islice(dicc.items(), 10))
# elementos = list(dicc.keys())
# conteos = list(dicc.values())
#
# df = pd.DataFrame()
# df['Nombre'] = elementos
# df['N veces'] = conteos
#
# df.to_excel('Bibliografía.xlsx')

for j in entrada:
    linia = j.rstrip().lstrip()
    if linia.startswith('<publisher>'):
        x = linia.replace('<publisher>', '')
        x = x.replace('</publisher>', '')
        lista.append(x)

counter = collections.Counter(lista)
dicc = dict(sorted(counter.items(), key=lambda item: item[1], reverse=True))
# diccionario = dict(itertools.islice(dicc.items(), 10))
elementos = list(dicc.keys())
conteos = list(dicc.values())

df = pd.DataFrame()
df['Nombre'] = elementos
df['N veces'] = conteos

df.to_excel('Bibliografía_editorial.xlsx')

# plt.figure(figsize=(50, 50))
# plt.bar(elementos, conteos)
# plt.title("Editorial", fontsize=50)
# plt.ylabel("Número de veces", fontsize=50)
# plt.xticks(rotation=90, ha="right", fontsize=10)
# plt.yticks(fontsize=50)
# plt.savefig('Bibliografía_editorial.jpg')
# plt.show()
#
