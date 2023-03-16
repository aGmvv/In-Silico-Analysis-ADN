import os
import tkinter
from itertools import product
from tkinter import messagebox
from tkinter import filedialog
from reverso_complementario import is_seq_valid


def create_combinations(secuencia):
    bases = []
    for base in secuencia:
        if base == "A":
            bases.append(["A"])
        elif base == "T":
            bases.append(["T"])
        elif base == "C":
            bases.append(["C"])
        elif base == "G":
            bases.append(["G"])
        elif base == "S":
            bases.append(["C", "G"])
        elif base == "W":
            bases.append(["A", "T"])
        elif base == "R":
            bases.append(["A", "G"])
        elif base == "Y":
            bases.append(["C", "T"])
        elif base == "K":
            bases.append(["G", "T"])
        elif base == "M":
            bases.append(["A", "C"])
        elif base == "B":
            bases.append(["C", "G", "T"])
        elif base == "D":
            bases.append(["A", "G", "T"])
        elif base == "H":
            bases.append(["A", "C", "T"])
        elif base == "V":
            bases.append(["A", "C", "G"])
        elif base == "N":
            bases.append(["A", "C", "G", "T"])

    # Generar todas las posibles combinaciones de bases
    combinaciones = list(product(*bases))

    # Convertir las combinaciones de bases en secuencias de ADN
    secuencias = ["".join(comb) for comb in combinaciones]

    save_file(secuencias)


def combinations_inicio(seq):
    seq = seq.upper()
    valor = is_seq_valid(seq)
    if valor:
        create_combinations(seq)
    else:
        print(f'La secuencia {seq} no es valida')


def save_file(secuencias):
    file_path = filedialog.asksaveasfilename(defaultextension=".txt",
                                             filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")])
    if file_path:
        if os.path.exists(file_path):
            response = tkinter.messagebox.askyesno("File already exists", "Do you want to overwrite the existing file?")
            if response == tkinter.YES:
                with open(file_path, "w") as file:
                    # Escribe cada elemento de la lista en una nueva línea del archivo
                    for elemento in secuencias:
                        file.write(elemento + '\n')
        else:
            with open(file_path, "w") as file:
                # Escribe cada elemento de la lista en una nueva línea del archivo
                for elemento in secuencias:
                    file.write(elemento + '\n')
