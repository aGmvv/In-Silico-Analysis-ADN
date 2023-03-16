"""
coding=utf-8
author = Alejandro González

Este script utiliza la biblioteca Selenium para automatizar el proceso de búsqueda de una taxonomía específica en la
base de datos del NCBI y la descarga de las referencias que coinciden con los criterios de búsqueda.

El script comienza iniciando un navegador Firefox y navegando a la URL especificada. A continuación, usa el método
find_element de Selenium para localizar e interactuar con varios elementos de la página web, como botones y campos de
texto, para realizar la búsqueda y navegar por los resultados.

El script también utiliza las bibliotecas Entrez y SeqIO del módulo Bio para interactuar con la base de datos NCBI y
recuperar la información de referencia. El script usa la biblioteca tqdm para mostrar una barra de progreso mientras
se ejecuta. Por último, el script guarda las referencias en un archivo especificado y cierra el navegador.
"""

import time
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio import Entrez
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options


def descargar_referencias(lista_ref, driver, number, por):
    """

    :param por: Porcentaje
    :param lista_ref: lista donde se van a añadir las Referencias que luego se descargaran.
    :param driver: el driver usado para acceder de forma remota a Internet.
    :param number: el tamaño de la primera secuencia, este tamaño se utiliza como referencia para calcular el 80%.

    Obtiene las secuencias y si la longitud es superior al 80% de la secuencia más larga la guarda en una lista.
    :return:
    """
    max_number = 1000000
    number = int(number * float(por))  # El 80% de la longitud de la secuencia más larga

    try:
        for i in range(max_number):
            try:
                problem = driver. find_element(By.XPATH, '/html/body/div[3]/div[2]/div/div[3]/button[2]')
                problem.click()

            except:
                # Obtiene el texto que se muestra en la web y lo separa para quedarse con las referencias.
                referencias = driver.find_element(By.XPATH, '/html/body/div[1]/div[1]/form/div[1]/div[4]/div/div[5]')
                # Separar el texto elemento a elemento es una lista usando '\n' como separador
                texto = referencias.text.split(sep='\n')

                for j in range(len(texto)):
                    # ------------------------------------------------------------------------------------------------ #
                    #                             Guardar la referencia si cumple las condiciones                      #
                    # ------------------------------------------------------------------------------------------------ #
                    if texto[j].startswith('Accession'):  # Comprobar si empieza por 'Accession'
                        try:
                            z = texto[j-1]
                            mayor = z.split(' ')  # Separar el string'' y obtener el primer elemento
                            # Convertir el primer elemento en un número y eliminar ','
                            mayor = int(mayor[0].replace(',', ''))
                            if mayor >= number:
                                # Añadir el elemento que empieza por 'Accession' y cumple con la distancia
                                # seleccionada a lista_ref
                                lista_ref.append(texto[j])
                            else:
                                return
                        except:
                            pass
                    # ----------------------------------------------------------------------------------------------- #

                time.sleep(2)

                # ---------------------------------------------------------------------------------------------------- #
                #                                         Pasar de página                                              #
                # ---------------------------------------------------------------------------------------------------- #
                if i == 0:  # Clic el boton next dependiendo de la iteración
                    next_button = driver.find_element(By.XPATH, "/html/body/div[1]/div[1]/form/div[1]/div[4]/div/div[3]"
                                                                "/div[2]/a[1]")
                    next_button.click()

                elif 0 < i < max_number:
                    next_button = driver.find_element(By.XPATH, "/html/body/div[1]/div[1]/form/div[1]/div[4]/div/div[3]"
                                                                "/div[2]/a[3]")
                    next_button.click()

                elif i == max_number:
                    break
                # ---------------------------------------------------------------------------------------------------- #
                time.sleep(30)
    except:
        return


def buscar(url, archivo, special_num, por):
    """

    :param por: Porcentaje
    :param url: enlace de NCBI.
    :param archivo: el nombre del archivo que se va a crear.
    :param special_num: número del taxid.

    Esta función abre un navegador de forma remota en el enlace que le hemos pasado en la variable url, filtra por
    taxonomy. Hay que tener en cuenta que no siempre el taxonomy que queremos es el primero, dependiendo del patogeno
    hay que cambiarlo.

    Después de seleccionar el taxonomy correcto, selecciona el número de secuencias que se muestran a 200 en la misma
    página y despues ordena las secuencias por orden de tamaño.

    Obtiene las referencias de todas las secuencias cuya longuitud es superior al 80% de la longuitud de la secuencia
    más larga.

    Después de obtener las referencias las descarga y guarda en un archivo fasta.
    :return:
    """
    options = Options()  # Añadir opciones al buscador
    options.add_argument("start-maximized")
    options.add_argument("disable-infobars")
    options.add_argument("--disable-extensions")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--no-sandbox")
    options.add_experimental_option('excludeSwitches', ['enable-logging'])

    driver = webdriver.Firefox()  # Crear un driver de Firefox
    driver.get(url)  # Abrir la URL
    tic = time.perf_counter()
    time.sleep(3)

    # ---------------------------------------------------------------------------------------------------------------- #
    #                                               Seleccionar taxonomy                                               #
    # ---------------------------------------------------------------------------------------------------------------- #
    taxonomy = driver.find_element(By.XPATH,
                                   f"/html/body/div[1]/div[1]/form/div[1]/div[5]/div/div[2]/div[2]/div/div/dl["
                                   f"2]/dd/div[1]/div[{special_num}]/span/a")
    taxonomy.click()
    time.sleep(10)

    # ---------------------------------------------------------------------------------------------------------------- #

    # ---------------------------------------------------------------------------------------------------------------- #
    #                                  Número de secuencias que se muestran por página                                 #
    # ---------------------------------------------------------------------------------------------------------------- #
    page = driver.find_element(By.XPATH, '/html/body/div[1]/div[1]/form/div[1]/div[4]/div/div[1]/ul/li[2]/a')
    page.click()
    time.sleep(1)

    number = driver.find_element(By.XPATH,
                                 '/html/body/div[1]/div[1]/form/div[1]/div[4]/div/div[1]/div[2]/fieldset/ul/li[6]/'
                                 'input')
    number.click()
    time.sleep(1)
    # ---------------------------------------------------------------------------------------------------------------- #
    # ---------------------------------------------------------------------------------------------------------------- #
    #                                        Ordenar por número de secuencia                                           #
    # ---------------------------------------------------------------------------------------------------------------- #
    sort = driver.find_element(By.XPATH, '/html/body/div[1]/div[1]/form/div[1]/div[4]/div/div[1]/ul/li[3]/a')
    sort.click()
    time.sleep(2)

    nsort = driver.find_element(By.XPATH,
                                '/html/body/div[1]/div[1]/form/div[1]/div[4]/div/div[1]/div[3]/fieldset/ul/li[7]/input')
    nsort.click()
    time.sleep(1)
    # ---------------------------------------------------------------------------------------------------------------- #
    tam = driver.find_element(By.XPATH, '/html/body/div[1]/div[1]/form/div[1]/div[4]/div/div[5]/div[1]/div[2]/div[1]/p')
    final_tam = tam.text.split(' ')
    number = final_tam[0]
    number = int(number.replace(',', ''))  # Obtener el tamaño de la secuencia más larga

    time.sleep(5)
    lista = []
    descargar_referencias(lista, driver, number, por)

    # ---------------------------------------------------------------------------------------------------------------- #
    #                  Separar las referencias para quedarnos solo con la información que nos interesa                 #
    # ---------------------------------------------------------------------------------------------------------------- #

    df = pd.DataFrame(lista)

    lista_ref = []
    df.columns = ['Referencias']

    for i in df['Referencias']:
        palabra = i.split(' ')
        lista_ref.append(palabra[1])
    # ---------------------------------------------------------------------------------------------------------------- #

    # ---------------------------------------------------------------------------------------------------------------- #
    #                      Descarga de NCBI las secuencias y las guarda en un arhivo .fasta                            #
    # ---------------------------------------------------------------------------------------------------------------- #
    Entrez.email = "agonzalez@certest.es"

    fasta = archivo + '.fasta'
    fw = open(fasta, 'w')

    for i, c in zip(lista_ref, tqdm(range(len(lista_ref)))):
        try:
            hd1 = Entrez.efetch(db="nucleotide", id=[i], rettype='fasta')
            seq = SeqIO.read(hd1, 'fasta')
            SeqIO.write(seq, fw, 'fasta')
        except:
            pass

    # Cierra el archivo una vez acabado
    fw.close()
    # ---------------------------------------------------------------------------------------------------------------- #

    print("FINISH")
    toc = time.perf_counter()
    print(f"Code in {toc - tic:0.4f} seconds")
