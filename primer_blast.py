import os
import time
import tkinter
from tqdm import tqdm
from Bio import SeqIO
from Bio import Entrez
from tkinter import filedialog
from tkinter import messagebox
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options


def buscar(url, records):
    options = Options()  # Añadir opciones al buscador
    options.add_argument("start-maximized")
    options.add_argument("disable-infobars")
    options.add_argument("--disable-extensions")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--no-sandbox")
    options.add_experimental_option('excludeSwitches', ['enable-logging'])

    driver = webdriver.Firefox()  # Crear un driver de Firefox
    driver.get(url)  # Abrir la URL

    time.sleep(3)

    first_primer = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[2]/div[1]/input')
    first_primer.send_keys(records[0].seq)

    time.sleep(1)

    second_primer = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/'
                                                  'fieldset[2]/div[2]/input')
    second_primer.send_keys(records[1].seq)

    time.sleep(1)

    database = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]/div['
                                             '3]/span/select')
    database.click()

    time.sleep(1)

    nr = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]/div['
                                       '3]/span/select/option[4]')
    nr.click()

    time.sleep(1)

    Organism = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]/div['
                                             '5]/div/div/input[1]')
    Organism.clear()

    time.sleep(1)

    pss1 = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]'
                                         '/div[7]/span[1]/select')
    pss1.click()

    time.sleep(1)

    pss11 = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]/div[7]/span['
                                          '1]/select/option[6]')
    pss11.click()

    time.sleep(1)

    pss2 = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]/'
                                         'div[7]/span[2]/select')
    pss2.click()

    time.sleep(1)

    pss12 = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]/div[7]/span['
                                          '2]/select/option[6]')
    pss12.click()

    time.sleep(1)

    amplicon = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]/div[8]/input')
    amplicon.clear()
    amplicon.send_keys('6000')

    time.sleep(1)

    get_primers = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/form/div[3]/div[1]/input')
    get_primers.click()

    time.sleep(1)
    valor = True

    while valor:
        try:
            texto = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[2]/div/div[3]/div')
            # Dividir el contenido del archivo en bloques separados por el carácter ">"
            bloques = texto.text.split(">")

            lista = []
            referencias = []
            inicio = []
            fin = []
            # Tiene que empezar a contar a partir del dos
            for i in range(2, len(bloques)):
                elemento = bloques[i].split('\n')
                lista.append((elemento[0], elemento[3], elemento[6]))

            for i in lista:
                for j in range(len(i)):
                    if j == 0:
                        referencias.append(i[j].split(' ')[0])
                    if j == 1:
                        inicio.append(i[j].split(' ')[8])
                    if j == 2:
                        fin.append(i[j].split(' ')[8])

            return referencias, inicio, fin

        except:
            time.sleep(10)


def primer_blast(archivo):
    tic = time.perf_counter()

    url = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/'
    records = list(SeqIO.parse(archivo, "fasta"))

    referencia, inicio_lista, fin_lista = buscar(url, records)

    # ---------------------------------------------------------------------------------------------------------------- #
    #                        Descarga de NCBI las secuencias y las guarda en un arhivo .fasta                          #
    # ---------------------------------------------------------------------------------------------------------------- #
    Entrez.email = "agonzalez@certest.es"

    fasta = save_file()
    fasta2 = save_file()

    fw = open(fasta, 'w')
    fw_acortamiento = open(fasta2, 'w')

    for i, inicio, fin, c in zip(referencia, inicio_lista, fin_lista, tqdm(range(len(referencia)))):
        try:
            hd1 = Entrez.efetch(db="nucleotide", id=[i], rettype='fasta')
            seq = SeqIO.read(hd1, 'fasta')
            SeqIO.write(seq, fw, 'fasta')
            if int(inicio) < int(fin):
                SeqIO.write(seq[int(inicio):int(fin)], fw_acortamiento, 'fasta')
            elif int(inicio) > int(fin):
                SeqIO.write(seq[int(fin):int(inicio)], fw_acortamiento, 'fasta')
        except:
            pass

    # Cierra el archivo una vez acabado
    fw.close()
    fw_acortamiento.close()

    print("FINISH")
    toc = time.perf_counter()
    print(f"Code in {toc - tic:0.4f} seconds")


def save_file():
    file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt"),
                                                                                 ("All Files", "*.*")])
    if file_path:
        if os.path.exists(file_path):
            response = tkinter.messagebox.askyesno("File already exists", "Do you want to overwrite the existing file?")
            if response == tkinter.YES:
                with open(file_path, "w") as file:
                    file.write("File content")
        else:
            with open(file_path, "w") as file:
                file.write("File content")

    return file_path
