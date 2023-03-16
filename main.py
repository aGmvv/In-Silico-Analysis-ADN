"""
coding=utf-8
author = Alejandro González

AGM
"""
import os
import sys
import Align
import tkinter
import warnings
import webbrowser
import circular_DNA
from LNA import LNA
from tkinter import *
from tkinter import ttk
from Blast import blast
from NCBI import buscar
from tabulate import tabulate
from tkinter import messagebox
from PIL import Image, ImageTk
from ID_primers import ID_primers
from find_nn_primers import contar
from Inclusivity import inclusivity
from tkinter import filedialog as fd
from primer_blast import primer_blast
from filtrar_sec import filtrar_archivo
from reverso_complementario import main_v2
from Cross_reactivity import preparar_cross
from primers_2_pr import primers_4_with_2_pr
from generate_combinations import combinations_inicio
from mutation_differentation import mutation_differentation


class Insilico:
    def __init__(self, master):
        # Crear variables
        self.url = None
        self.bar = None
        self.name = None
        self.extrar = None
        self.extrap = None
        self.evalue = None
        self.maxseqs = None
        self.nodiana = None
        self.master = master
        self.numdiana = None
        self.taxidnum = None
        self.pregunta = None
        self.etiquetas = None
        self.botonstep = None
        self.botonerror = None
        self.entryerror = None
        self.etiqueta31 = None
        self.etiquetand = None
        self.botonstep2 = None
        self.etiquetaen = None
        self.etiquetafc = None
        self.etiquetana = None
        self.porcentaje = None
        self.casilla_out = None
        self.prueba_name = None
        self.player_name = None
        self.etiquetaflo = None
        self.etiquetauno = None
        self.NCBI_ventana = None
        self.etiqueta_pat = None
        self.player_name3 = None
        self.player_name4 = None
        self.etiquetana33 = None
        self.etiquetana32 = None
        self.player_name1 = None
        self.player_name2 = None
        self.etiquetana34 = None
        self.otra_ventana = None
        self.ventana_error = None
        self.etiquetaerror = None
        self.blast_ventana = None
        self.etiquetablast = None
        self.otra_ventanalo = None
        self.etiquetarerror = None
        self.player_name_pat = None
        self.casilla_reverso = None
        self.etiquetareverso = None
        self.reverso_ventana = None
        self.circular_ventana = None
        self.etiquetataxidnum = None
        self.botonfuncionblast = None
        self.etiquetaporcentraje = None
        self.pregunta_inclusivity = None

        # --------------------------------------------------------------------------------------------------------- #
        #                                            Configurar menu                                                #
        # --------------------------------------------------------------------------------------------------------- #
        menu = Menu(self.master)
        filemenu = Menu(menu)
        filemenu.add_command(label='Count mutations', command=self.count_mutations)
        filemenu.add_command(label='get references', command=get_references)
        menu.add_cascade(label="Tools", menu=filemenu)

        align_menu = Menu(menu)
        align_menu.add_command(label='Aling', command=funcion_align)
        menu.add_cascade(label='Aling', menu=align_menu)

        special_menu = Menu(menu)
        special_menu.add_command(label='Circular DNA - primers', command=circular_dna)
        special_menu.add_command(label='4 primers (2pr)', command=self.primers_4)
        menu.add_cascade(label='Special case', menu=special_menu)

        ID = Menu(menu)
        ID.add_command(label='Reverse', command=self.reverso)
        ID.add_command(label='See IUPAC', command=see_IUPAC)
        ID.add_command(label='Primer Blast', command=primerblast)
        ID.add_command(label='Primer list', command=ID_primers)
        ID.add_command(label='Combinations ADN', command=self.create_combinations)
        ID.add_command(label='Most popular nn', command=most_popular)
        menu.add_cascade(label='I+D', menu=ID)

        option_menu = Menu(menu)
        option_menu.add_command(label='Exit', command=root.destroy)
        menu.add_cascade(label='Options', menu=option_menu)

        self.master.config(menu=menu)
        # --------------------------------------------------------------------------------------------------------- #
        """ 
        
        Pruebas es un comando que asigna archivos para hacer pruebas, los archivos estan en la carpeta de datos 
        de prueba. Esta función la he creado para no tener que meter todo el rato los datos a mano cuando quiero ver
        que uno de los cambios que estoy haciendo funciona o no.
        
        Generalmente esta siempre comentada para que no aparezca en la aplicación, solo la des comento cuando voy 
        a hacer cambios.
        
        Esta funcion SOLO esta para el analysis inclusivity.
        """
        # filemenu.add_command(label='Pruebas', command=self.prueba_inclusivity)

        # --------------------------------------------------------------------------------------------------------- #
        #                                         Pestaña de inicio                                                 #
        # --------------------------------------------------------------------------------------------------------- #
        """                            Los tres botones (Inclusivity, Blast y Cross)                              """

        if hasattr(sys, '_MEIPASS'):
            path = sys._MEIPASS
        else:
            path = os.getcwd()

        image1 = Image.open(os.path.join(path, "Sources/Logo_nuevo.png"))
        image1 = image1.resize((150, 50))
        self.test = ImageTk.PhotoImage(image1)
        self.web = Button(master, command=abrir_pagina_web, image=self.test, compound=TOP, bg='white')
        self.web.place(relx=0.01, rely=0.01)

        self.etiqueta = Label(master, text="IN SILICO ANALYSIS", font=("Century Gothic", 30, "bold"), bg='white')
        self.etiqueta.pack()

        self.nombre = Label(master, text="AGM", font=("Century Gothic", 8), bg='white')
        self.nombre.place(relx=0.97, rely=0.97)

        image2 = Image.open(os.path.join(path, "Sources/inclusivity.jpg"))
        image2 = image2.resize((200, 100))
        self.imagen_inclusivity = ImageTk.PhotoImage(image2)
        self.botoncroos = Button(master, text="Inclusivity", command=self.inclusivity,
                                 image=self.imagen_inclusivity, compound=TOP, bg='#c6ced1')
        self.botoncroos.place(relx=0.2, rely=0.1)

        image4 = Image.open(os.path.join(path, "Sources/blast.jpg"))
        image4 = image4.resize((200, 100))
        self.imagen_blast = ImageTk.PhotoImage(image4)
        self.botonblast = Button(master, text="Blast", command=self.funcion_blast,
                                 image=self.imagen_blast, compound=TOP, bg='#c6ced1')
        self.botonblast.place(relx=0.45, rely=0.1)

        image3 = Image.open(os.path.join(path, "Sources/cross.jpg"))
        image3 = image3.resize((200, 100))
        self.imagen_cross = ImageTk.PhotoImage(image3)
        self.botonlootro = Button(master, text="Cross Reactivity", command=self.funcion_cross,
                                  image=self.imagen_cross, compound=TOP, bg='#c6ced1')
        self.botonlootro.place(relx=0.7, rely=0.1)

        self.mmboton = Button(master, text="Alineamientos", command=self.funcion_inclusivity, bg='#c6ced1')
        self.mmboton.place(relx=0.587, rely=0.226)

        self.separador = ttk.Separator(master, orient='horizontal')
        self.separador.pack(fill='x', padx=1, pady=250)

        # --------------------------------------------------------------------------------------------------------- #

    def count_mutations(self):
        """

        :return:
        """
        self.otra_ventana = tkinter.Toplevel(root)
        self.otra_ventana.wm_title("In silico Analysis")
        self.otra_ventana.geometry("500x400")

        self.etiquetafc = Label(self.otra_ventana, text="Count mutations")
        self.etiquetafc.pack()

        # Casilla
        self.etiquetaen = Label(self.otra_ventana, text="Nombre del archivo que se va a generar")
        self.etiquetaen.pack()

        self.player_name = ttk.Entry(self.otra_ventana)
        self.player_name.insert(END, ".xlsx")
        self.player_name.get()
        self.player_name.pack()

        self.botonstep = Button(self.otra_ventana, text="Enter", command=self.step_count_mutations)
        self.botonstep.place(x=230, y=100)

    def step_count_mutations(self):
        name = self.player_name.get()
        count_mutations(name)

        self.etiquetas = Label(self.otra_ventana, text="El archivo ha sido generado")
        self.etiquetas.pack()

    def primers_4(self):
        """

        Esta función corresponde en el menu de Special case cuando queremos alinear un conjunto de 4 primers teniendo
        dos sondas.

        Fw - Pr1 - Pr2 - Rv
        :return:
        """
        self.otra_ventana = tkinter.Toplevel(root)
        self.otra_ventana.wm_title("In silico Analysis")
        self.otra_ventana.geometry("500x400")

        self.etiquetafc = Label(self.otra_ventana, text="4 primers, 2 sondas")
        self.etiquetafc.pack()

        # Casilla
        self.etiquetaen = Label(self.otra_ventana, text="Nombre del archivo que se va a generar")
        self.etiquetaen.pack()

        self.player_name = ttk.Entry(self.otra_ventana)
        self.player_name.insert(END, ".xlsx")
        self.player_name.get()
        self.player_name.pack()

        self.botonstep = Button(self.otra_ventana, text="Enter", command=self.step_primers_4)
        self.botonstep.place(x=230, y=100)

    def step_primers_4(self):
        name = self.player_name.get()
        datos_primers_4(name)

        self.etiquetas = Label(self.otra_ventana, text="El archivo ha sido generado")
        self.etiquetas.pack()

    def inclusivity(self):
        """

        La función corresponde al boton inclusivity, a esta función hay que proporcionarle el enlace de NCBI del
        patogeno que queremos buscar.
        Por ejemplo: https://www.ncbi.nlm.nih.gov/nuccore/?term=zika

        Además, hay que añadir el nombre del archivo que se va a crear, como se indican en las etiquetas.
        :return:
        """
        self.NCBI_ventana = tkinter.Toplevel(root)
        self.NCBI_ventana.wm_title("Inclusivity")
        self.NCBI_ventana.geometry("500x400")

        self.etiquetablast = Label(self.NCBI_ventana, text="Introduce el enlace de NCBI")
        self.etiquetablast.pack()

        self.url = ttk.Entry(self.NCBI_ventana)
        self.url.insert(END, "")
        self.url.get()
        self.url.pack()

        self.etiquetablast = Label(self.NCBI_ventana, text="El nombre del archivo que se va a crear")
        self.etiquetablast.pack()

        self.name = ttk.Entry(self.NCBI_ventana)
        self.name.insert(END, "")
        self.name.get()
        self.name.pack()

        self.etiquetataxidnum = Label(self.NCBI_ventana, text="Posición del patógeno en la columna by taxon")
        self.etiquetataxidnum.pack()

        self.taxidnum = ttk.Entry(self.NCBI_ventana)
        self.taxidnum.insert(END, "1")
        self.taxidnum.get()
        self.taxidnum.pack()

        self.etiquetaporcentraje = Label(self.NCBI_ventana, text="Porcentaje")
        self.etiquetaporcentraje.pack()

        self.porcentaje = ttk.Entry(self.NCBI_ventana)
        self.porcentaje.insert(END, "0.8")
        self.porcentaje.get()
        self.porcentaje.pack()

        self.botonfuncionblast = Button(self.NCBI_ventana, text="Enter", command=self.final_NCBI)
        self.botonfuncionblast.place(x=230, y=205)

    def final_NCBI(self):
        """

        Guarda en variables el enlace añadido y el nombre del archivo que se va a generar.
        :return:
        """
        url = self.url.get()
        archivo = self.name.get()
        special_num = self.taxidnum.get()
        por = self.porcentaje.get()
        buscar(url, archivo, special_num, por)

    def funcion_blast(self):
        """

        Esta función corresponde al boton Blast. Es un paso intermedio para seleccionar el archivo que se va a crear.
        :return:
        """
        self.blast_ventana = tkinter.Toplevel(root)
        self.blast_ventana.wm_title("Blast")
        self.blast_ventana.geometry("500x400")

        self.etiquetablast = Label(self.blast_ventana, text="Introduce el nombre del archivo que se va a crear")
        self.etiquetablast.pack()

        self.casilla_out = ttk.Entry(self.blast_ventana)
        self.casilla_out.insert(END, ".txt")
        self.casilla_out.get()
        self.casilla_out.pack()

        self.etiquetablast = Label(self.blast_ventana, text="e-value")
        self.etiquetablast.pack()

        self.evalue = ttk.Entry(self.blast_ventana)
        self.evalue.insert(END, "10000")
        self.evalue.get()
        self.evalue.pack()

        self.etiquetablast = Label(self.blast_ventana, text="Número de secuencias máximas")
        self.etiquetablast.pack()

        self.maxseqs = ttk.Entry(self.blast_ventana)
        self.maxseqs.insert(END, "8000")
        self.maxseqs.get()
        self.maxseqs.pack()

        self.etiquetablast = Label(self.blast_ventana, text="Puntuación por match")
        self.etiquetablast.pack()

        self.extrar = ttk.Entry(self.blast_ventana)
        self.extrar.insert(END, "1")
        self.extrar.get()
        self.extrar.pack()

        self.etiquetablast = Label(self.blast_ventana, text="Penalización por missmatch")
        self.etiquetablast.pack()

        self.extrap = ttk.Entry(self.blast_ventana)
        self.extrap.insert(END, "-1")
        self.extrap.get()
        self.extrap.pack()

        self.botonfuncionblast = Button(self.blast_ventana, text="Enter", command=self.final_blast)
        self.botonfuncionblast.place(x=230, y=240)

    def final_blast(self):
        """

        Esta función se divide en dos partes.
        1.- Introducir los primers que vamos a usar para hacer el blast.
        2.- El archivo que se genera en el blast es el .txt que va a utilizar para realizar el siguiente paso del Cross.
            Cuando se acaba el paso del blast, en la pestaña creada con Tkinter se muestra un boton que deja de
            continuar automáticamente con el siguiente paso, aunque de momento, hay que volver a seleccionar los
            archivos.
        :return:
        """
        tkinter.messagebox.showinfo('title', 'Introducir los primers')
        primers = fd.askopenfilename()

        patogeno = self.casilla_out.get()
        evalue = self.evalue.get()
        max_seqs = self.maxseqs.get()
        extrar = self.extrar.get()
        extrap = self.extrap.get()

        blast(primers, patogeno, evalue, max_seqs, extrar, extrap)

        self.pregunta_inclusivity = Button(self.blast_ventana, text="Continuar con el primer paso de Cross",
                                           command=self.funcion_inclusivity)
        self.pregunta_inclusivity.place(x=100, y=210)

    def create_combinations(self):
        """

        Esta función realiza el reverso complementario de la secuencia que se le proporciona. Tenemos en cuenta las
        letras especiales que corresponden a los nucleotidos degenerados.

        Esta función es un paso intermedio en la que se introduce la secuencia de ADN.
        :return:
        """
        self.reverso_ventana = tkinter.Toplevel(root)
        self.reverso_ventana.wm_title("Generate combinations of ADN sequences")
        self.reverso_ventana.geometry("500x400")

        self.etiquetareverso = Label(self.reverso_ventana, text="Introduce la cadena de ADN")
        self.etiquetareverso.pack()

        self.casilla_reverso = ttk.Entry(self.reverso_ventana)
        self.casilla_reverso.insert(END, "")
        self.casilla_reverso.get()
        self.casilla_reverso.pack()

        self.botoncroos = Button(self.reverso_ventana, text="Enter", command=self.main_combinations)
        self.botoncroos.place(x=230, y=70)

    def main_combinations(self):
        """

        Guarda la secuencia de ADN y la envia al archivo reverso_complementario.py, donde se calcula el reverso
        complementario.
        :return:
        """
        sec = self.casilla_reverso.get()
        combinations_inicio(sec)

    def reverso(self):
        """

        Esta función realiza el reverso complementario de la secuencia que se le proporciona. Tenemos en cuenta las
        letras especiales que corresponden a los nucleotidos degenerados.

        Esta función es un paso intermedio en la que se introduce la secuencia de ADN.
        :return:
        """
        self.reverso_ventana = tkinter.Toplevel(root)
        self.reverso_ventana.wm_title("Reverse Complement")
        self.reverso_ventana.geometry("500x400")

        self.etiquetareverso = Label(self.reverso_ventana, text="Introduce la cadena de ADN")
        self.etiquetareverso.pack()

        self.casilla_reverso = ttk.Entry(self.reverso_ventana)
        self.casilla_reverso.insert(END, "")
        self.casilla_reverso.get()
        self.casilla_reverso.pack()

        self.botoncroos = Button(self.reverso_ventana, text="Enter", command=self.main_reverso)
        self.botoncroos.place(x=230, y=70)

    def main_reverso(self):
        """

        Guarda la secuencia de ADN y la envia al archivo reverso_complementario.py, donde se calcula el reverso
        complementario.
        :return:
        """
        sec = self.casilla_reverso.get()
        main_v2(sec)

    def step_inclusivity(self):
        """

        Aunque la función se llame step_inclusivity, es el paso intermedio del boton alineamiento.
        Cuando se empezó el código esto correspondía al paso de inclusivity, cuando se fue desarrollando el código y
        este cambio, para no tener que cambiar las funciones y tener problemas de código se dejó así.
        :return:
        """
        name = self.player_name.get()
        datos_inclusivity(name)

        self.etiquetas = Label(self.otra_ventana, text="El archivo ha sido generado")
        self.etiquetas.pack()

    def paso_intermedio(self):
        """

        Paso intermedio para el paso cross, corresponde al boton cross en la aplicación. En esta función se selecciona
        el número de dianas.
        :return:
        """
        diana = self.numdiana.get()
        diana = int(diana)

        if diana == 1:
            name1 = self.player_name1.get()
            datos_cross1(name1, diana)
        if diana == 2:
            name1 = self.player_name1.get()
            name2 = self.player_name2.get()
            datos_cross2(name1, name2, diana)
        if diana == 3:
            name1 = self.player_name1.get()
            name2 = self.player_name2.get()
            name3 = self.player_name3.get()
            datos_cross3(name1, name2, name3, diana)
        if diana == 4:
            name1 = self.player_name1.get()
            name2 = self.player_name2.get()
            name3 = self.player_name3.get()
            name4 = self.player_name4.get()
            datos_cross4(name1, name2, name3, name4, diana)

        self.etiqueta = Label(self.otra_ventanalo, text="Los archivos han sido generados")
        self.etiqueta.pack()

    def dianas(self):
        """

        Selección de archivos dependiendo el número de dianas. Paso Cross.
        :return:
        """
        diana = self.numdiana.get()
        diana = int(diana)
        if diana == 1:
            # Casilla
            self.etiquetauno = Label(self.otra_ventanalo, text="Nombre del primer archivo")
            self.etiquetauno.pack()

            self.player_name1 = ttk.Entry(self.otra_ventanalo)
            self.player_name1.insert(END, ".xlsx")
            self.player_name1.get()
            self.player_name1.pack()

            self.botonstep = Button(self.otra_ventanalo, text="Enter", command=self.paso_intermedio)
            self.botonstep.pack()

        if diana == 2:
            # Primera casilla
            self.etiqueta = Label(self.otra_ventanalo, text="Nombre del primer archivo")
            self.etiqueta.pack()

            self.player_name1 = ttk.Entry(self.otra_ventanalo)
            self.player_name1.insert(END, ".xlsx")
            self.player_name1.get()
            self.player_name1.pack()

            # Segunda casilla
            self.etiquetana = Label(self.otra_ventanalo, text="Nombre del segundo archivo")
            self.etiquetana.pack()

            self.player_name2 = ttk.Entry(self.otra_ventanalo)
            self.player_name2.insert(END, ".xlsx")
            self.player_name2.get()
            self.player_name2.pack()

            self.botonstep = Button(self.otra_ventanalo, text="Enter", command=self.paso_intermedio)
            self.botonstep.pack()

        if diana == 3:
            # Primera casilla
            self.etiqueta31 = Label(self.otra_ventanalo, text="Nombre del primer archivo")
            self.etiqueta31.pack()

            self.player_name1 = ttk.Entry(self.otra_ventanalo)
            self.player_name1.insert(END, ".xlsx")
            self.player_name1.get()
            self.player_name1.pack()

            # Segunda casilla
            self.etiquetana32 = Label(self.otra_ventanalo, text="Nombre del segundo archivo")
            self.etiquetana32.pack()

            self.player_name2 = ttk.Entry(self.otra_ventanalo)
            self.player_name2.insert(END, ".xlsx")
            self.player_name2.get()
            self.player_name2.pack()

            # Tercera casilla
            self.etiquetana33 = Label(self.otra_ventanalo, text="Nombre del tercer archivo")
            self.etiquetana33.pack()

            self.player_name3 = ttk.Entry(self.otra_ventanalo)
            self.player_name3.insert(END, ".xlsx")
            self.player_name3.get()
            self.player_name3.pack()

            self.botonstep = Button(self.otra_ventanalo, text="Enter", command=self.paso_intermedio)
            self.botonstep.pack()

        if diana == 4:
            # Primera casilla
            self.etiqueta31 = Label(self.otra_ventanalo, text="Nombre del primer archivo")
            self.etiqueta31.pack()

            self.player_name1 = ttk.Entry(self.otra_ventanalo)
            self.player_name1.insert(END, ".xlsx")
            self.player_name1.get()
            self.player_name1.pack()

            # Segunda casilla
            self.etiquetana32 = Label(self.otra_ventanalo, text="Nombre del segundo archivo")
            self.etiquetana32.pack()

            self.player_name2 = ttk.Entry(self.otra_ventanalo)
            self.player_name2.insert(END, ".xlsx")
            self.player_name2.get()
            self.player_name2.pack()

            # Tercera casilla
            self.etiquetana33 = Label(self.otra_ventanalo, text="Nombre del tercer archivo")
            self.etiquetana33.pack()

            self.player_name3 = ttk.Entry(self.otra_ventanalo)
            self.player_name3.insert(END, ".xlsx")
            self.player_name3.get()
            self.player_name3.pack()

            # Cuarta casilla
            self.etiquetana34 = Label(self.otra_ventanalo, text="Nombre del cuarto archivo")
            self.etiquetana34.pack()

            self.player_name4 = ttk.Entry(self.otra_ventanalo)
            self.player_name4.insert(END, ".xlsx")
            self.player_name4.get()
            self.player_name4.pack()

            self.botonstep = Button(self.otra_ventanalo, text="Enter", command=self.paso_intermedio)
            self.botonstep.pack()

    def funcion_inclusivity(self):
        """

        Esta función corresponde al boton alineamiento debajo del boton de Blast.
        :return:
        """
        self.otra_ventana = tkinter.Toplevel(root)
        self.otra_ventana.wm_title("In silico Analysis")
        self.otra_ventana.geometry("500x400")

        self.etiquetafc = Label(self.otra_ventana, text="Primer paso cross")
        self.etiquetafc.pack()

        # Casilla
        self.etiquetaen = Label(self.otra_ventana, text="Nombre del archivo que se va a generar")
        self.etiquetaen.pack()

        self.player_name = ttk.Entry(self.otra_ventana)
        self.player_name.insert(END, ".xlsx")
        self.player_name.get()
        self.player_name.pack()

        self.botonstep = Button(self.otra_ventana, text="Enter", command=self.step_inclusivity)
        self.botonstep.place(x=230, y=70)

    def funcion_cross(self):
        """

        Pestaña desplegable correspondiente al boton Cross de la aplicación.
        :return:
        """
        self.otra_ventanalo = tkinter.Toplevel(root)
        self.otra_ventanalo.wm_title("In silico Analysis")
        self.otra_ventanalo.geometry("500x400")

        self.etiquetaflo = Label(self.otra_ventanalo, text="Cross Reactivity")
        self.etiquetaflo.pack()

        self.etiquetand = Label(self.otra_ventanalo, text="Seleccionar número de dianas, entre 1 y 4")
        self.etiquetand.pack()

        self.numdiana = ttk.Entry(self.otra_ventanalo)
        self.numdiana.insert(END, "1")
        self.numdiana.get()
        self.numdiana.pack()

        self.botonstep = Button(self.otra_ventanalo, text="Enter", command=self.dianas)
        self.botonstep.place(x=230, y=70)


def most_popular():
    tkinter.messagebox.showinfo('title', 'Introducir archivo')
    oligos = fd.askopenfilename()
    contar(oligos)


def primerblast():
    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers')
    oligos = fd.askopenfilename()  # Introducir los primers en formato fasta
    primer_blast(oligos)


def prueba_inclusivity():
    """

    Esta función es de prueba, se usa para no tener que estar seleccionando los archivos una y otra vez cuando se
    prueba una modificación en el código.

    :return:
    """
    excel_name = ""
    filename = open("")
    oligos = open("")

    # Llamar a la función inclusivity para probarla con los archivos seleccionados.
    inclusivity(filename, oligos, excel_name)


def get_references():
    tkinter.messagebox.showinfo('title', 'Introducir el archivo a analizar')
    filename = fd.askopenfilename()
    filtrar_archivo(filename)


def see_IUPAC():
    """

    :return:
    """
    headers = ["", "A", "C", "G", "T"]
    table = [["Y", "-", "C", "-", "T"],
             ["W", "A", "-", "-", "T"],
             ["S", "-", "C", "G", "-"],
             ["M", "A", "C", "-", "-"],
             ["K", "-", "-", "G", "T"],
             ["R", "A", "-", "G", "-"],
             ["B", "-", "C", "G", "T"],
             ["D", "A", "-", "G", "T"],
             ["H", "A", "C", "-", "T"],
             ["V", "A", "C", "G", "-"],
             ["N", "A", "C", "G", "T"]]

    print(tabulate(table, headers, tablefmt="rounded_outline"))


def circular_dna():
    """

    Esta función corresponde a la pestaña del menu Special case, Circular DNA - primers.
    El archivo Align.py y circular_DNA.py son iguales, lo unico que cambia entre los dos archivos es que en
    circular_DNA.py se duplica la secuencia para "resolver" el problema de alinear los primers con secuencias de ADN
    circular.
    :return:
    """
    circular_DNA.inicio()


def funcion_align():
    """

    Esta función llama al archivo Align.py.
    :return:
    """
    Align.inicio()


def datos_LNA(excel_name, patogeno):
    """

    :param excel_name: nombre del archivo que se va a generar.
    :param patogeno: nombre del patogeno que queremos buscar.

    Esta función está en desarrollo y actualmente NO está en uso.
    :return:
    """
    filename = 'C:/Users/agonzalez/Desktop/Tkinter_v2.py/In_silico/lna/listeria.txt'
    oligos = 'C:/Users/agonzalez/Desktop/Tkinter_v2.py/In_silico/lna/listeria.fasta'

    LNA(filename, oligos, excel_name, patogeno)


def count_mutations(excel_name):
    tkinter.messagebox.showinfo('title', 'Introducir el archivo a analizar')
    filename = fd.askopenfilename()  # El primer archivo es el .txt que obtenemos al hacer el blast

    mutation_differentation(filename, excel_name)


def datos_primers_4(excel_name):
    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo a analizar')
    filename = fd.askopenfilename()  # El primer archivo es el .txt que obtenemos al hacer el blast

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers')
    oligos = fd.askopenfilename()  # Introducir los primers en formato fasta

    primers_4_with_2_pr(filename, oligos, excel_name)


def datos_inclusivity(excel_name):
    """

    :param excel_name: nombre del archivo que se va a generar

    Esta función corresponde al boton alineamiento, se encarga de seleccionar los archivos que vamos a usar.
    :return:
    """
    # Seleccionar archivos para Inclusivity
    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo a analizar')
    filename = fd.askopenfilename()  # El primer archivo es el .txt que obtenemos al hacer el blast

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers')
    oligos = fd.askopenfilename()  # Introducir los primers en formato fasta

    inclusivity(filename, oligos, excel_name)


def datos_cross1(name1, diana):
    """

    :param name1: nombre del primer archivo que se va a generar, en ese caso el unico.
    :param diana: número de diana.

    Cuando seleccionas el número de diana 1, introduces el nombre del archivo que se va a analizar y los primers, como
    el resto de dianas están vacías se completan como nulas.
    :return:
    """
    # Seleccionar archivos para lo otro
    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo a analizar')
    filename = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 1')
    oligos1 = fd.askopenfilename()

    name2, name3, name4 = '', '', ''
    oligos2, oligos3, oligos4 = '', '', ''

    preparar_cross(filename, name1, name2, name3, name4, oligos1, oligos2, oligos3, oligos4, diana)


def datos_cross2(name1, name2, diana):
    """

    :param name1: nombre del primer archivo que se va a generar.
    :param name2: nombre del segundo archivo que se va a generar.
    :param diana: número de diana.

    Cuando seleccionas el número de diana 2, introduces el nombre de los archivos que se va a analizar y los primers,
    como el resto de dianas están vacías se completan como nulas.
    :return:
    """
    # Seleccionar archivos para lo otro
    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo a analizar')
    filename = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 1')
    oligos1 = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 2')
    oligos2 = fd.askopenfilename()

    name3, name4 = '', ''
    oligos3, oligos4 = '', ''

    preparar_cross(filename, name1, name2, name3, name4, oligos1, oligos2, oligos3, oligos4, diana)


def datos_cross3(name1, name2, name3, diana):
    """

    :param name1: nombre del primer archivo que se va a generar.
    :param name2: nombre del segundo archivo que se va a generar.
    :param name3: nombre del tercer archivo que se va a generar.
    :param diana: número de diana.

    Cuando seleccionas el número de diana 3, introduces el nombre de los archivos que se va a analizar y los primers,
    como el resto de dianas están vacías se completan como nulas.
    :return:
    """
    # Seleccionar archivos para lo otro
    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo a analizar')
    filename = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 1')
    oligos1 = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 2')
    oligos2 = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 3')
    oligos3 = fd.askopenfilename()

    name4 = ''
    oligos4 = ''

    preparar_cross(filename, name1, name2, name3, name4, oligos1, oligos2, oligos3, oligos4, diana)


def datos_cross4(name1, name2, name3, name4, diana):
    """

    :param name1: nombre del primer archivo que se va a generar.
    :param name2: nombre del segundo archivo que se va a generar.
    :param name3: nombre del tercer archivo que se va a generar.
    :param name4: nombre del cuarto archivo que se va a generar.
    :param diana: número de diana.

    Cuando seleccionas el número de diana 3, introduces el nombre de los archivos que se va a analizar y los primers,
    como el resto de dianas están vacías se completan como nulas.
    :return:
    """
    # Seleccionar archivos para lo otro
    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo a analizar')
    filename = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 1')
    oligos1 = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 2')
    oligos2 = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 3')
    oligos3 = fd.askopenfilename()

    tkinter.messagebox.showinfo('title', 'Introducir nombre del archivo con los primers 4')
    oligos4 = fd.askopenfilename()

    preparar_cross(filename, name1, name2, name3, name4, oligos1, oligos2, oligos3, oligos4, diana)


def abrir_pagina_web():
    webbrowser.open_new("https://www.certest.es/")


if __name__ == '__main__':
    warnings.filterwarnings('ignore')  # Evitar que se muestren los warning

    # ------------------------------------------------------------#
    #        Inicio de la aplicación y llamada a la función       #
    # ------------------------------------------------------------#
    root = Tk()
    miVentana = Insilico(root)
    root.wm_title("In Silico Analysis")
    root.configure(bg='white')
    root.state('zoomed')

    root.config()
    root.mainloop()
    # ------------------------------------------------------------#
