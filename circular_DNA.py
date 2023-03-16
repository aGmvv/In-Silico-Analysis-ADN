import os
import csv
import time
import ali_ui
import operator
import itertools
import xlsxwriter
import pandas as pd
from itertools import zip_longest
from reverso_complementario import main
from PyQt5.QtWidgets import QApplication, QComboBox, QMainWindow, QTabWidget, QFileDialog


def send_mail(to, file):
    from email.mime.multipart import MIMEMultipart
    from email.mime.text import MIMEText
    import smtplib
    # create message object instance
    msg = MIMEMultipart()

    message = "Alignment " + file + " finished"

    # setup the parameters of the message
    password = "4L1gn3r2021"
    msg['From'] = "aligner@certest.bio"
    msg['To'] = to
    msg['Subject'] = "Align:"

    # add in the message body
    msg.attach(MIMEText(message, 'plain'))

    # create server
    server = smtplib.SMTP_SSL('smtp.serviciodecorreo.es', 465)

    # Login Credentials for sending the mail
    server.login(msg['From'], password)
    try:
        # send the message via the server
        server.sendmail(msg['From'], msg['To'], msg.as_string())
        print("successfully sent email to %s:" % (msg['To']))
    except:
        print("Not E-Mail specified")
    finally:
        server.quit()


def order_two_list(list_a, list_b):
    temp = list(zip(list_a, list_b))
    temp.sort(key=operator.itemgetter(0))
    list_a, list_b = zip(*temp)
    return list_a, list_b


def enumerate_files(path):
    # --------------------------------------------------------------------------------------------- #
    #                                       Guardar el direcotrio                                   #
    # --------------------------------------------------------------------------------------------- #
    print(path)
    dirname = path
    files = os.listdir(dirname)
    flist = []
    for f in files:
        file = path + '/' + f
        flist.append(file)
    print("FILO:", flist)
    # --------------------------------------------------------------------------------------------- #
    return flist


def calc(OlAlign, v):
    score = OlAlign.score
    path = OlAlign[v].path
    query = OlAlign[v].query[
            path[0][1] - path[0][0]:path[1][1] + (len(OlAlign[0].target) - path[1][0])]
    return score, path, query


def csv_writer(data):
    with open("align.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_wr = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_wr.writerows(data)


def compare_seq(target, query, oli_name, multiple, bnc, one_and_wrong, pr):
    # --------------------------------------------------------------------------------------------- #
    #                                  Dar los valores numericos                                    #
    # --------------------------------------------------------------------------------------------- #
    result = ''
    positions = []
    counter = 0
    valor = 0
    for count, (t, q) in enumerate(zip(target, query)):
        if t == q or q == 'N':  # Si es igual o es una N se asigna como homologo
            result = result + '-'
            counter = counter + 1
        else:
            if oli_name == 'Forward' or oli_name == 'Reverse':  # Solo para forward o reverse
                if multiple == 'No':  # Saber si es multiple
                    if count < 5 or count >= len(query) - 5:  # Si el mm esta entre los nn de 0 a 5 en cualquier extremo
                        valor = 3  # Codigo 3; grave
                    else:
                        valor = 2  # Codigo 2; se detecta pero hay retroceso en la grafica
                else:  # Si es multiple se pone 4
                    valor = 4
            elif oli_name == 'Probe':  # La sonda
                if pr in one_and_wrong:  # Si la sonda esta mutada en alguno de estos kits ya esta mal
                    valor = 3  # Codigo 3; grave
                elif multiple == 'No':
                    if 5 < count <= len(query) - 5:  # Si el mm esta entre los nn de 0 a 5 en cualquier extremo
                        valor = 3  # Codigo 3; grave
                    else:
                        valor = 2  # Codigo 2; se detecta pero hay retroceso en la grafica
                else:  # Si es multiple se pone 4
                    valor = 4
            if q in bnc:  # bnc son las letras especiales
                if t in bnc[q]:  # Si hay una de estas letras tambien se detecta como homologo
                    result = result + '-'
                    counter = counter + 1
                else:
                    result = result + q
            elif t in bnc:
                if q in bnc[t]:
                    result = result + '-'
                    counter = counter + 1
                else:
                    result = result + q
            else:
                result = result + q

            positions.append(counter)
    mm = (len(target) - counter)
    if mm >= 3 and multiple == 'No':
        valor = 3
    return result, mm, valor


def excel(workbook, res, chains, total, target1, target2):
    worksheet = workbook.add_worksheet('Resume')
    worksheet3 = workbook.add_worksheet('Resume2')
    col = 1
    row = 2
    merge_format = workbook.add_format({'bold': 1, 'border': 1, 'align': 'center', 'valign': 'vcenter'})

    # Worksheet1
    worksheet.merge_range('B1:D1', target1[5], merge_format)
    worksheet.merge_range('E1:G1', target2[5], merge_format)
    worksheet.write('A2', 'Record')
    worksheet.write('B2', 'FW')
    worksheet.write('C2', 'RV')
    worksheet.write('D2', 'SD')
    worksheet.write('E2', 'FW2')
    worksheet.write('F2', 'RV2')
    worksheet.write('G2', 'SD2')
    worksheet.write('H2', 'Total')
    for at, reg in enumerate(chains):
        for iteration, cha in enumerate(reg):
            worksheet.write(row, col, cha)
            worksheet.write(row, 7, total[at])
            col = col + 1
        col = 1
        row = row + 1
    # Worksheet3
    worksheet3.write('A1', 'Both 100% Homology')
    worksheet3.write('B1', 'Target 1 100% Homology')
    worksheet3.write('C1', 'Target 2 100% Homology')
    worksheet3.write('D1', 'Both <100% Homology')
    worksheet3.write('E1', 'Target 1 <100% Homology')
    worksheet3.write('F1', 'Target 2 <100% Homology')
    row = 1
    col = 0
    for tot in res:
        worksheet3.write(row, col, tot)
        col = col + 1


def excel2(worksheet2, score, mm, path, result, target, record, variant, oliname, pr, query, control, registro,
           multiple, assays_name):
    column_name = ["-", "-", "-", "-"]
    for i, name in enumerate(assays_name):
        column_name[i] = name

    if control == 0:
        worksheet2.write('A1', 'Registro')
        worksheet2.write('B1', 'Oliname')
        worksheet2.write('C1', 'Ensayo')
        worksheet2.write('D1', 'Score')
        worksheet2.write('E1', 'MM')
        worksheet2.write('F1', 'Inicio')
        worksheet2.write('G1', 'Fin')
        worksheet2.write('H1', 'Target')
        worksheet2.write('I1', 'Query')
        worksheet2.write('J1', 'Result')
        worksheet2.write('K1', 'Variant')
        worksheet2.write('L1', 'Chain')
        worksheet2.write('M1', '¿Multiple?')
        worksheet2.write('N1', 'CHAIN2')
        worksheet2.write('O1', 'SeqGen')
        worksheet2.write('P1', 'MultiClass')
        worksheet2.write('Q1', 'ComboVariant')
        worksheet2.write('R1', 'SingleClass')
        worksheet2.write('S1', 'Max ' + column_name[0])
        worksheet2.write('T1', 'Max ' + column_name[1])
        worksheet2.write('U1', 'Max ' + column_name[2])
        worksheet2.write('V1', 'Max ' + column_name[3])
        control = 1

    worksheet2.write(registro, 0, record)
    worksheet2.write(registro, 1, oliname)
    worksheet2.write(registro, 2, pr)
    worksheet2.write(registro, 3, score)
    worksheet2.write(registro, 4, mm)
    worksheet2.write(registro, 5, path[0][1])
    worksheet2.write(registro, 6, path[1][1])
    worksheet2.write(registro, 7, str(target))
    worksheet2.write(registro, 8, str(query))
    worksheet2.write(registro, 9, result)
    worksheet2.write(registro, 10, variant)
    worksheet2.write(registro, 12, multiple)
    registro = registro + 1
    return control, registro


def write_chain(registro, registro_anterior, cadena, worksheet2):
    listToStr = ''.join([str(elem) for elem in cadena])
    for x in range(registro_anterior, registro):
        worksheet2.write(x, 11, listToStr)


def write_chain2(registro, registro_anterior, cadena, worksheet2, len_align, variant, assays_name, df_corr):
    i = 1
    listToStr = ''.join([str(elem) for elem in cadena])

    # Comprobamos si esta NCO5 y NCO11
    if "NCO5" and "NCO11" in assays_name:
        cl, cl2, reason = cover(listToStr, assays_name, df_corr)
    else:
        cl, cl2 = classifier(listToStr)
    for x in range(registro_anterior, registro):
        worksheet2.write(x, 13, listToStr)
        worksheet2.write(x, 15, cl)
        worksheet2.write(x, 16, variant[:-1])
        worksheet2.write(x, 17, cl2)
        if len(len_align) == 1:
            if i in range(1, len_align[0] + 1):
                worksheet2.write(x, 14, listToStr[0:3])
                worksheet2.write(x, 18, max(cadena[0:3]))
                i = i + 1
        if len(len_align) == 2:

            if i in range(1, len_align[0] + 1):
                worksheet2.write(x, 14, listToStr[0:3])
                worksheet2.write(x, 18, max(cadena[0:3]))
                worksheet2.write(x, 19, max(cadena[3:6]))
            if i in range(len_align[0] + 1, len_align[0] + len_align[1] + 1):
                worksheet2.write(x, 14, listToStr[3:6])
                worksheet2.write(x, 18, max(cadena[0:3]))
                worksheet2.write(x, 19, max(cadena[3:6]))
            i = i + 1
        if len(len_align) == 3:
            if i in range(1, len_align[0] + 1):
                worksheet2.write(x, 14, listToStr[0:3])
                worksheet2.write(x, 18, max(cadena[0:3]))
                worksheet2.write(x, 19, max(cadena[3:6]))
                worksheet2.write(x, 20, max(cadena[6:9]))

            if i in range(len_align[0] + 1, len_align[0] + len_align[1] + 1):
                worksheet2.write(x, 14, listToStr[3:6])
                worksheet2.write(x, 18, max(cadena[0:3]))
                worksheet2.write(x, 19, max(cadena[3:6]))
                worksheet2.write(x, 20, max(cadena[6:9]))

            if i in range(len_align[0] + len_align[1] + 1, len_align[0] + len_align[1] + len_align[2] + 1):
                worksheet2.write(x, 14, listToStr[6:9])
                worksheet2.write(x, 18, max(cadena[0:3]))
                worksheet2.write(x, 19, max(cadena[3:6]))
                worksheet2.write(x, 20, max(cadena[6:9]))
            i = i + 1
        if len(len_align) == 4:
            if i in range(1, len_align[0] + 1):
                worksheet2.write(x, 14, listToStr[0:3])
                worksheet2.write(x, 18, max(cadena[0:3]))
                worksheet2.write(x, 19, max(cadena[3:6]))
                worksheet2.write(x, 20, max(cadena[6:9]))
                worksheet2.write(x, 21, max(cadena[9:12]))

            if i in range(len_align[0] + 1, len_align[0] + len_align[1] + 1):
                worksheet2.write(x, 14, listToStr[3:6])
                worksheet2.write(x, 18, max(cadena[0:3]))
                worksheet2.write(x, 19, max(cadena[3:6]))
                worksheet2.write(x, 20, max(cadena[6:9]))
                worksheet2.write(x, 21, max(cadena[9:12]))

            if i in range(len_align[0] + len_align[1] + 1, len_align[0] + len_align[1] + len_align[2] + 1):
                worksheet2.write(x, 14, listToStr[6:9])
                worksheet2.write(x, 18, max(cadena[0:3]))
                worksheet2.write(x, 19, max(cadena[3:6]))
                worksheet2.write(x, 20, max(cadena[6:9]))
                worksheet2.write(x, 21, max(cadena[9:12]))

            if i in range(len_align[0] + len_align[1] + len_align[2] + 1,
                          len_align[0] + len_align[1] + len_align[2] + len_align[3] + 1):
                worksheet2.write(x, 14, listToStr[9:12])
                worksheet2.write(x, 18, max(cadena[0:3]))
                worksheet2.write(x, 19, max(cadena[3:6]))
                worksheet2.write(x, 20, max(cadena[6:9]))
                worksheet2.write(x, 21, max(cadena[9:12]))
            i = i + 1


def classifier(string):
    strlen = len(string)
    max_compare = 0
    max_compare_single = 0

    if strlen == 3:
        first3 = list(string[0:3])
        if int(max(first3)) != 0:
            max_compare = int(max(first3))
        max_compare_single = int(max(first3))
    if strlen == 6:
        first3 = list(string[0:3])
        second3 = list(string[3:6])
        if int(max(first3)) != 0 and int(max(second3)) != 0:
            max_compare = min(int(max(first3)), int(max(second3)))
        max_compare_single = min(int(max(first3)), int(max(second3)))
    if strlen == 9:
        first3 = list(string[0:3])
        second3 = list(string[3:6])
        third3 = list(string[6:9])

        if int(max(first3)) != 0 and int(max(second3)) != 0 and int(max(third3)) != 0:
            # print('MAX1= ', max(first3), ' MAX2= ', max(second3))
            max_compare = min(int(max(first3)), int(max(second3)), int(max(third3)))
        max_compare_single = min(int(max(first3)), int(max(second3)), int(max(third3)))
    if strlen == 12:
        first3 = list(string[0:3])
        second3 = list(string[3:6])
        third3 = list(string[6:9])
        fourth3 = list(string[9:12])

        if int(max(first3)) != 0 and int(max(second3)) != 0 and int(max(third3)) != 0 and int(max(fourth3)) != 0:
            max_compare = min(int(max(first3)), int(max(second3)), int(max(third3)), int(max(fourth3)))
        max_compare_single = min(int(max(first3)), int(max(second3)), int(max(third3)), int(max(fourth3)))

    return max_compare, max_compare_single


def chain_cut(cadena, n, fillvalue='x'):
    args = [iter(cadena)] * n
    ans = list(zip_longest(fillvalue=fillvalue, *args))
    t = len(ans)
    for i in range(t):
        ans[i] = "".join(ans[i])
    return ans


def cover(cadena, assays_name, df_corr):
    mmax = 333
    assays_number = len(assays_name)
    indexes = {}
    total_list = []
    l = []
    a = 1
    c_cut = chain_cut(cadena, assays_number)
    for i in range(0, len(cadena)):
        l.append(cadena[i])
        if a % 3 == 0:
            total_list.append(l)
            l = []
        a += 1
    for nb, jk in zip(assays_name, total_list):
        indexes[nb] = jk

    nco5index = assays_name.index("NCO5")
    nco11index = assays_name.index("NCO11")

    nco5loc = str(c_cut[nco5index])
    nco11loc = str(c_cut[nco11index])
    othernco = max(total_list[0])

    locs = df_corr.loc[(df_corr['NCO5'] == nco5loc) & (df_corr['NCO11'] == nco11loc)]

    if othernco == '0':
        mmax = locs.loc[:, 0].values[0]
    if othernco == '2':
        mmax = locs.loc[:, 2].values[0]
    if othernco == '3':
        mmax = locs.loc[:, 3].values[0]
    if othernco == '4':
        mmax = locs.loc[:, 4].values[0]

    return mmax, mmax, ''


def write_chain_csv(list):
    list_to_str = ''.join([str(elem) for elem in list])
    return list_to_str


def store_value(valor, chain):
    c = chain.append(valor)
    return c


def mutation(target1, target2, target3, target4, seq, oliname):
    variant = ''
    targets = [target1, target2, target3, target4]
    for tg in targets:
        if tg[1] == 'NCO1':
            if oliname == 'Reverse':
                if seq == 'TCC--GAG-------------':
                    variant = 'Cluster V'
        if tg[1] == 'NCO2' and tg[5] == 'NCO5':
            if oliname == 'Forward':
                if seq == '------T---------------':
                    variant = 'Sudafrica'
                if seq == 'AAC-------------------':
                    variant = 'AAC Mutation'
                if seq == 'T---------------------':
                    variant = 'India'
            if oliname == 'Reverse':
                if seq == '-----------------T----':
                    variant = 'Brasil'
                if seq == '-------------------T--':
                    variant = 'UK-3T'
        if tg[1] == 'NCO2' and tg[5] == 'NCO9':
            if oliname == 'Probe':
                if seq == '------------A----------':
                    variant = 'P13A'
        if tg[5] == 'N501Y':
            if oliname == 'Forward':
                if seq == '-----------------A----':
                    variant = 'N501Y_F_1'
        if tg[5] == 'NCO7':
            if oliname == 'Probe':
                if seq == '--T---------------------':
                    variant = 'omicron'
                if seq == '--T------------------G--':
                    variant = 'BA.5'
        if tg[5] == 'NCO9':
            if oliname == 'Forward':
                if seq == '---------------A----':
                    variant = 'BA.4'
        if tg[5] == 'NCO4':
            if oliname == 'Probe':
                if seq == '----T-----------------------':
                    variant = 'P.1.2'
            if oliname == 'Forward':
                if seq == 'A-AATG---------------':
                    variant = 'B.1.1.138'
    return variant


def parameters(ndiana, target1, target2, target3, target4, ford, reve, sond, ford2, reve2, sond2, ford3, reve3,
               sond3, ford4, reve4, sond4):
    """ --- Decimos el numero de dianas y le pasamos los valores --- """
    global chain3
    oli = []
    oli_name = []
    assay = 1
    assays_name = []
    if ndiana == "1":
        assay = 1
        assays_name = [target1[5]]
        oli_name = ["Forward", "Reverse", "Probe"]
        oli = [ford, reve, sond]
        chain3 = [0, 0, 0]
    if ndiana == "2":
        assay = 2
        assays_name = [target1[5], target2[5]]
        oli_name = ["Forward", "Reverse", "Probe", "Forward", "Reverse", "Probe"]
        oli = [ford, reve, sond, ford2, reve2, sond2]
        chain3 = [0, 0, 0, 0, 0, 0]
    if ndiana == "3":
        assay = 3
        assays_name = [target1[5], target2[5], target3[5]]
        oli_name = ["Forward", "Reverse", "Probe", "Forward", "Reverse", "Probe", "Forward", "Reverse", "Probe"]
        oli = [ford, reve, sond, ford2, reve2, sond2, ford3, reve3, sond3]
        chain3 = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    if ndiana == "4":
        assay = 4
        assays_name = [target1[5], target2[5], target3[5], target4[5]]
        oli_name = ["Forward", "Reverse", "Probe", "Forward", "Reverse", "Probe", "Forward", "Reverse", "Probe",
                    "Forward", "Reverse", "Probe"]
        oli = [ford, reve, sond, ford2, reve2, sond2, ford3, reve3, sond3, ford4, reve4, sond4]
        chain3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    return assay, assays_name, oli_name, oli, chain3


class MainWindow(QMainWindow, QComboBox, QApplication, QTabWidget, ali_ui.Ui_MainWindow):

    def __init__(self):
        """ La función super() accede a los métodos heredados anulados en una clase. La función super() se usa
        en la clase hija con herencia múltiple para acceder a la función de la siguiente clase padre o superclase.
        La función super() determina la siguiente clase principal utilizando el Orden de resolución del método (MRO).
        Como si el MRO es C -> D -> B -> A -> object, para D, la función super() buscará la siguiente clase padre o
        método de superclase en la secuencia D -> B -> A -> object. """
        super(MainWindow, self).__init__()

        # ----------------------------------------------------------------------------------------------- #
        #                                  Crear entorno de la aplicacion                                 #
        # ----------------------------------------------------------------------------------------------- #
        self.filename = None
        self.ui = ali_ui.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.pushButtonFile.clicked.connect(self.browse_slot)
        self.ui.pushButtonFile_2.clicked.connect(self.browse_slot2)
        # self.ui.pushButtonFile_3.clicked.connect(self.xmlparser)
        self.ui.pushButton_2.clicked.connect(self.resume)
        # self.ui.pushButton_3.clicked.connect(self.blast)
        # ----------------------------------------------------------------------------------------------- #

        # ----------------------------------------------------------------------------------------------- #
        #                                   Iniciar primers y dianas                                      #
        # ----------------------------------------------------------------------------------------------- #
        self.ui.geneCombo.currentTextChanged.connect(lambda x: self.set_primers(x, "1"))
        self.ui.geneCombo_2.currentTextChanged.connect(lambda x: self.set_primers(x, "2"))
        self.ui.geneCombo_3.currentTextChanged.connect(lambda x: self.set_primers(x, "3"))
        self.ui.geneCombo_4.currentTextChanged.connect(lambda x: self.set_primers(x, "4"))
        self.ui.comboDianas.currentTextChanged.connect(self.setdianas)
        self.ui.lineMScore.setText("1")
        self.ui.lineMMScore.setText("-1")
        self.ui.lineOGScore.setText("-5")
        self.ui.lineEGScore.setText("-2")
        self.ui.linePos1.setText("0")
        self.ui.linePos2.setText("0")
        self.ui.linePos3.setText("0")
        self.ui.linePos4.setText("0")
        self.ui.comboBox.setCurrentIndex(6)
        # ----------------------------------------------------------------------------------------------- #

    def setdianas(self, s):
        # ----------------------------------------------------------------------------------------------- #
        #                       Seleccionar el numero de dianas que vamos a utlizar                       #
        # ----------------------------------------------------------------------------------------------- #
        fields2 = [self.ui.lineFW2, self.ui.lineRV2, self.ui.lineSD2, self.ui.linePos2, self.ui.geneCombo_2]
        fields3 = [self.ui.lineFW3, self.ui.lineRV3, self.ui.lineSD3, self.ui.linePos3, self.ui.geneCombo_3]
        fields4 = [self.ui.lineFW4, self.ui.lineRV4, self.ui.lineSD4, self.ui.linePos4, self.ui.geneCombo_4]
        if s == "1":
            for f in [*fields2, *fields3, *fields4]:
                f.setDisabled(True)
        if s == "2":
            for f in fields2:
                f.setEnabled(True)
            for f in [*fields3, *fields4]:
                f.setDisabled(True)
        if s == "3":
            for f in [*fields2, *fields3]:
                f.setEnabled(True)
            for f in [*fields4]:
                f.setDisabled(True)
        if s == "4":
            for f in [*fields2, *fields3, *fields4]:
                f.setEnabled(True)
        # ----------------------------------------------------------------------------------------------- #

    def set_primers(self, s, name):
        """
        Las parametros dados son los primers y el número del ensayo.
        Declarar los primers, s es el nombre de la diana. Los datos son los primers (Fw, Rv, Pr), añade información
        adicional.
        """
        # --------------------------------------------------------------------------------------------- #
        #          Seleccionamos los primes que queremos, en el caso de que los metamos a mano,         #
        #          seleccionamos la opcion                                                              #
        #          E1, ... y los metemos a mano                                                         #
        # --------------------------------------------------------------------------------------------- #
        data = []
        emptys = ['E1', 'E2', 'E3', 'E4']
        if s in emptys:
            data = ["", "", "", "", "", "", ""]

        if s == "NCO2":
            data = ["CCCCTGCATACACTAATTCTT", "AGACATGTATAGCATGGAACC", "CACACGTGGTGTTTATTACCCTGAC", "S", "NCO1", "21563",
                    "25384"]
        if s == "NCO2.5":
            data = ["GTGGTGTTTATTACCCTGAC", "GGGTTATCAAACCTCTTAGTAC", "TATACATGTCTCTGG", "S", "SUK1", "21563",
                    "25384"]
        if s == "NCO2.6":
            data = ["GTGGTGTTTATTACCCTGAC", "GGGTTATCAAACCTCTTAGTAC", "CCATGCTATCTCTGGG", "S", "SUK1/SUK2", "21563",
                    "25384"]
        if s == "NCO28":
            data = ["ACTTCTTTTTCTTGCTTTCGTG", "AGCAGTACGCACACAATC", "CTAGTTACACTAGCCATCCTTACTGC", "E", "NCO4", "26245",
                    "26472"]
        if s == "NCO4":
            data = ["CCCTGTGGGTTTTACACTTAA", "ACGATTGTGCATCAGCTGA", "CCGTCTGCGGTATGTGGAAAGGTTATGG", "ORF1ab", "NCO2",
                    "13000", "14000"]
        if s == "NCO5":
            data = ["GGGGAACTTCTCCTGCTAGAAT", "CAGACATTTTGCTCTCAAGCTG", "TTGCTGCTGCTTGACAGATT", "N", "NCO2",
                    "28274", "29533"]
        if s == "NCO7":
            data = ["GACCCCAAAATCAGCGAAAT", "TCTGGTTACTGCCAGTTGAATCTG", "ACCCCGCATTACGTTTGGTGGACC", "N", "NCO3",
                    "28274", "29533"]
        if s == "NCO9":
            data = ["TTACAAACATTGGCCGCAAA", "GCGCGACATTCCGAAGAA", "ACAATTTGCCCCCAGCGCTTCAG", "N", "NCO2",
                    "28274", "29533"]
        if s == "E484K":
            data = ["TCAAACCTTTTGAGAGAGATATT", "TGGAAACCATATGATTGTAAAGG", "ATGGTGTTAAAGG", "S", "VAR",
                    "21563", "25384"]
        if s == "E484Q":
            data = ["TCAAACCTTTTGAGAGAGATATT", "TGGAAACCATATGATTGTAAAGG", "TGTTCAAGGTTT", "S", "VAR",
                    "21563", "25384"]
        if s == "K417N":
            data = ["AGTCAGACAAATCGCTCCAG", "CGCAGCCTGTAAAATCATCT", "TGGAAATATTGCTGA", "S", "VAR",
                    "21563", "25384"]
        if s == "K417T":
            data = ["AGTCAGACAAATCGCTCCAG", "CGCAGCCTGTAAAATCATCT", "TGGAACGATTGC", "S", "VAR",
                    "21563", "25384"]
        if s == "L452R":
            data = ["CTTGATTCTAAGGTTGGTGGT", "CACCATTACAAGGTGTGCT", "TACCGGTATAGA", "S", "VAR",
                    "21563", "25384"]
        if s == "N501Y":
            data = ["CACCTTGTAATGGTGTTGAAGG", "ACTACTCTGTATGGTTGGTAAC", "CCCACTTATGGT", "S", "VAR",
                    "21563", "25384"]
        if s == "P681R":
            data = ["GTGCAGGTATATGCGCTAGTT", "AGTGTAGGCAATGATGGATTGA", "ATTCTCGTCGGC", "S", "VAR",
                    "21563", "25384"]
        if s == "NCO11":
            data = ["AAATTTTGGGGACCAGGAAC", "TGGCACCTGTGTAGGTCAAC", "ATGTCGCGCATTGGCATGGA", "N", "NCO2-Reforzado",
                    "28274", "29533"]
        if s == "YIA1":
            data = ["GACCRATCCTGTCACCTCTGAC", "AGGGCATTYTGGACAAAKCGTCTA", "CGTGCCCAGTGAGCGAGGACTGCA",
                    "Matrix Gene (M1)", "Influenza A",
                    "28274", "29533"]
        if s == "YIB1":
            data = ["GAGACACAATTGCCTACCTGCTT", "TTCTTTCCCACCGAACCAAC", "AGAAGATGGAGAAGGCAAAGCAGAACTAGC",
                    "Matrix Gene (M1)", "Influenza B",
                    "28274", "29533"]
        if s == "HNV1":
            data = ["GGCTGCTTTGAATTTTACCACAA", "TTTGGGTAGTCATAAGTCCCATTTT", "TGCGATAACACGTGCATGGAAAGTGTC",
                    "Hemaglutinina", "Influenza A, H1N1v",
                    "28274", "29533"]
        if s == "RSB1":
            data = ["GATGGCTCTTAGCAAAGTCAAGTTAA", "TGTCAATATTATCTCCTGTACTACGTTGAA",
                    "TGATACATTAAATAAGGATCAGCTGCTGTCATCCA", "N", "Respiratoty syncitial virus B",
                    "28274", "29533"]
        if s == "RSA1":
            data = ["GCTCTTAGCAAAGTCAAGTTGAATGA", "TGCTCCGTTGGATGGTGTATT", "ACACTCAACAAAGATCAACTTCTGTCATCCAGC", "N",
                    "Respiratoty syncitial virus A",
                    "28274", "29533"]
        if s == "VAR11":
            data = ["GTGACCTTGGTGCTTGTAT", "CGTAGTTGTTCAGACAATGAC", "CAACATTACTTTGA", "N", "VAO",
                    "28274", "29533"]
        if s == "VAR13":
            data = ["TTTAATAGTGCTATTGGCAAAATTC", "TTGGAGCTAAGTTGTTTAACAAG", "AACCATAATGCA", "N", "VAO",
                    "28274", "29533"]
        if s == "HSV1.2":
            data = ["TATTGGTGCGATGGCGACAC", "CTTTCCGCATGTGGGCTCTC", "GCGGGTAGGGTATGGGGCGGGG", "US4", "HHT",
                    "", ""]
        if s == "BDT.2":
            data = ["TCGAACGCGTGGAATGG", "GGCCGTTGGCTTCAAATAGA", "AGACCCAGGGCGCACGCTGTC", "pIS1001", "BDT",
                    "", ""]
        if s == "BDT.3":
            data = ["GGCGACAGCGAGACAGAATC", "GCCGCCTTGGCTCACTT", "CGTGCAGATAGGCTTTTAGCTTGAGCGC", "hIS1001", "BDT",
                    "", ""]

        """ 
        Llama a la funcion self.assign_primers y le pasa como parametros los primers y la informacion adicional y el 
        numero del ensayo
        """
        self.assign_primers(data, name)

    def assign_primers(self, data, number):
        # --------------------------------------------------------------------------------------------- #
        #                                       Guardar los primers                                     #
        # --------------------------------------------------------------------------------------------- #
        if number == "1":
            self.ui.lineEdit.setText(data[0])
            self.ui.lineEdit_2.setText(data[1])
            self.ui.lineEdit_3.setText(data[2])
            self.ui.lineGen1.setText(data[3])
            self.ui.lineProd1.setText(data[4])
            self.ui.lineInicio1.setText(data[5])
            self.ui.lineFin1.setText(data[6])
        if number == "2":
            self.ui.lineFW2.setText(data[0])
            self.ui.lineRV2.setText(data[1])
            self.ui.lineSD2.setText(data[2])
            self.ui.lineGen2.setText(data[3])
            self.ui.lineProd2.setText(data[4])
            self.ui.lineInicio2.setText(data[5])
            self.ui.lineFin2.setText(data[6])
        if number == "3":
            self.ui.lineFW3.setText(data[0])
            self.ui.lineRV3.setText(data[1])
            self.ui.lineSD3.setText(data[2])
            self.ui.lineGen3.setText(data[3])
            self.ui.lineProd3.setText(data[4])
            self.ui.lineInicio3.setText(data[5])
            self.ui.lineFin3.setText(data[6])
        if number == "4":
            self.ui.lineFW4.setText(data[0])
            self.ui.lineRV4.setText(data[1])
            self.ui.lineSD4.setText(data[2])
            self.ui.lineGen4.setText(data[3])
            self.ui.lineProd4.setText(data[4])
            self.ui.lineInicio4.setText(data[5])
            self.ui.lineFin4.setText(data[6])

    def browse_slot(self):
        # --------------------------------------------------------------------------------------------- #
        #                                   Seleccionar el archivo .fasta                               #
        # --------------------------------------------------------------------------------------------- #
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name = QFileDialog.getExistingDirectory()
        if file_name:
            print(file_name)
            self.ui.lineFile.setText(file_name)
        # --------------------------------------------------------------------------------------------- #
        return file_name

    def browse_slot2(self):
        """ Crear torno de la aplicacion """
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(None, "QFileDialog.getOpenFileName()", "",
                                                   "All Files (*);;Python Files (*.py)", options=options)
        if file_name:
            print(file_name)
            self.ui.lineFile_2.setText(file_name)
        return file_name

    def resume(self):
        from Bio.Seq import Seq

        tic = time.perf_counter()
        df_corr = pd.read_excel('C:/Users/agonzalez/Desktop/Tkinter_v2.py/In_silico/Sources/correla2.xlsx',
                                dtype=str, engine='openpyxl')

        # --------------------------------------------------------------------------------------------- #
        #                                  Darle la vuelta al reverse                                   #
        # --------------------------------------------------------------------------------------------- #
        rvseq = Seq(self.ui.lineEdit_2.text())
        rvseq2 = Seq(self.ui.lineRV2.text())
        rvseq3 = Seq(self.ui.lineRV3.text())
        rvseq4 = Seq(self.ui.lineRV4.text())
        # --------------------------------------------------------------------------------------------- #
        rvseqrvc = rvseq.reverse_complement()
        rvseqrvc2 = rvseq2.reverse_complement()
        rvseqrvc3 = rvseq3.reverse_complement()
        rvseqrvc4 = rvseq4.reverse_complement()
        ndiana = self.ui.comboDianas.currentText()
        # --------------------------------------------------------------------------------------------- #

        # --------------------------------------------------------------------------------------------- #
        #        Los target contienen la informacion adicional que hemos guardado al principio          #
        # --------------------------------------------------------------------------------------------- #
        target1 = [self.ui.lineGen1.text(), self.ui.lineProd1.text(), self.ui.lineInicio1.text(),
                   self.ui.lineFin1.text(), self.ui.linePos1.text(), self.ui.geneCombo.currentText()]

        target2 = [self.ui.lineGen2.text(), self.ui.lineProd2.text(), self.ui.lineInicio2.text(),
                   self.ui.lineFin2.text(), self.ui.linePos2.text(), self.ui.geneCombo_2.currentText()]

        target3 = [self.ui.lineGen3.text(), self.ui.lineProd3.text(), self.ui.lineInicio3.text(),
                   self.ui.lineFin3.text(), self.ui.linePos3.text(), self.ui.geneCombo_3.currentText()]

        target4 = [self.ui.lineGen4.text(), self.ui.lineProd4.text(), self.ui.lineInicio4.text(),
                   self.ui.lineFin4.text(), self.ui.linePos4.text(), self.ui.geneCombo_4.currentText()]
        # --------------------------------------------------------------------------------------------- #

        # --------------------------------------------------------------------------------------------- #
        #                                      Hacer el alineamiento                                    #
        # --------------------------------------------------------------------------------------------- #
        self.ali(self.ui.lineEdit.text(), rvseqrvc, self.ui.lineEdit_3.text(), self.ui.lineFW2.text(),
                 rvseqrvc2, self.ui.lineSD2.text(), self.ui.lineFW3.text(), rvseqrvc3,
                 self.ui.lineSD3.text(), self.ui.lineFW4.text(), rvseqrvc4, self.ui.lineSD4.text(),
                 self.ui.lineFile.text(), target1, target2, target3, target4, ndiana, df_corr)
        # --------------------------------------------------------------------------------------------- #
        # var_check = filtros_final(self.filename)
        #
        # var_check.to_excel(self.filename, sheet_name='Bad')

        print("FINISH")
        send_mail(self.ui.lineEmail.text(), self.ui.lineFile.text())
        toc = time.perf_counter()
        print(f"Code in {toc - tic:0.4f} seconds")

    def ali(self, ford, reve, sond, ford2, reve2, sond2, ford3, reve3, sond3, ford4, reve4, sond4, file2, target1,
            target2, target3, target4, ndiana, df_corr):
        from Bio import SeqIO
        from Bio import Align

        one_and_wrong = ['E484K', 'N501Y', 'K417T', 'K417N', 'NCO2.6', 'E484Q', 'P681R', 'L452R', 'VAR11', 'VAR13']

        pr = ''
        file = enumerate_files(file2)

        # --------------------------------------------------------------------------------------------- #
        #                                      Empieza el alineamiento                                  #
        # --------------------------------------------------------------------------------------------- #
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = int(self.ui.lineMScore.text())
        aligner.mismatch_score = int(self.ui.lineMMScore.text())
        aligner.open_gap_score = int(self.ui.lineOGScore.text())
        aligner.extend_gap_score = int(self.ui.lineEGScore.text())
        # --------------------------------------------------------------------------------------------- #

        for fi in file:
            records = list(SeqIO.parse(fi, "fasta"))  # Lee el archivo fasta
            self.filename = fi.replace('.fasta', '') + ".xlsx"  # Crea el nombre del archivo .xlsx
            bnc = {'Y': ['C', 'T'],
                   'W': ['A', 'T'],
                   'S': ['C', 'G'],
                   'M': ['A', 'C'],
                   'K': ['G', 'T'],
                   'R': ['A', 'G'],
                   'B': ['C', 'G', 'T'],
                   'D': ['A', 'G', 'T'],
                   'H': ['A', 'C', 'T'],
                   'V': ['A', 'C', 'G']}
            workbook = xlsxwriter.Workbook(self.filename)
            worksheet2 = workbook.add_worksheet('Bad')  # Crea el archivo y la pestaña Bad

            assay, assays_name, oli_name, oli, chain3 = parameters(ndiana, target1, target2, target3, target4,
                                                                   ford, reve, sond, ford2, reve2, sond2, ford3,
                                                                   reve3, sond3, ford4, reve4, sond4)

            control = 0
            registro = 1

            targs1 = [0, 1, 2]
            targs2 = [3, 4, 5]
            targs3 = [6, 7, 8]
            targs4 = [9, 10, 11]
            registro_anterior = 1

            mm = 0
            tmm = []
            num = 2

            for i in range(0, len(records)):
                variant_chain = ''
                cadena = []
                lengths = []

                for d in range(0, int(ndiana)):
                    lengths.append(0)

                print(i + 1, "-", records[i].name)  # Aparece por pantalla el numero y el nombre de la secuencia

                for a, x in enumerate(oli):
                    if a in targs1:
                        pr = target1[5]

                    if a in targs2:
                        pr = target2[5]

                    if a in targs3:
                        pr = target3[5]

                    if a in targs4:
                        pr = target4[5]

                    lista_secuencias = list(itertools.repeat(str(records[i].seq), num))
                    cycle = ''.join(lista_secuencias)

                    # --------------------------------------------------------------------------------------------- #
                    #                              Alinea y asigna numero de gravedad                               #
                    # --------------------------------------------------------------------------------------------- #
                    OlAlign = aligner.align(x, cycle.upper())
                    OlAlign_reverso = aligner.align(main(x), cycle.upper())

                    """ Comprobar si el primer esta en sentido contrario a la secuencia """
                    if OlAlign_reverso.score > OlAlign.score:
                        OlAlign = OlAlign_reverso

                    if len(OlAlign) == 1:  # Si es 1 es porque no tiene alineamientos multiples
                        n = 0
                        multiple = "No"
                        if len(OlAlign[n].target) == OlAlign[n].score:  # Quiere decir que alinea al 100%
                            tmm.append(0)
                            cadena.append(0)
                        else:
                            cadena.append(1)
                            tmm.append(mm)
                    else:  # Si no es 1 tiene alineamientos multiples
                        cadena.append(2)
                        multiple = "Yes"

                    for n in range(len(OlAlign)):
                        score, path, query = calc(OlAlign, n)

                        result, mm, valor = compare_seq(OlAlign[n].target, query, oli_name[a], multiple, bnc,
                                                        one_and_wrong, pr)

                        variant = mutation(target1, target2, target3, target4, result, oli_name[a])
                        if variant != '':  # Si la mutacion esta experimentalmente comprobada
                            valor = 1
                            variant_chain = variant_chain + variant + '-'

                        chain3[a] = valor
                        if mm == 0:
                            valor = 0
                        chain3[a] = valor

                        if len(lengths) == 1:
                            if pr == assays_name[0]:
                                lengths[0] = lengths[0] + 1
                        if len(lengths) == 2:
                            if pr == assays_name[0]:
                                lengths[0] = lengths[0] + 1
                            if pr == assays_name[1]:
                                lengths[1] = lengths[1] + 1

                        if len(lengths) == 3:
                            if pr == assays_name[0]:
                                lengths[0] = lengths[0] + 1
                            if pr == assays_name[1]:
                                lengths[1] = lengths[1] + 1
                            if pr == assays_name[2]:
                                lengths[2] = lengths[2] + 1

                        if len(lengths) == 4:
                            if pr == assays_name[0]:
                                lengths[0] = lengths[0] + 1
                            if pr == assays_name[1]:
                                lengths[1] = lengths[1] + 1
                            if pr == assays_name[2]:
                                lengths[2] = lengths[2] + 1
                            if pr == assays_name[3]:
                                lengths[3] = lengths[3] + 1

                        control, registro = excel2(worksheet2, score, mm, path, result, OlAlign[n].target,
                                                   records[i].name, variant, oli_name[a], pr, query, control,
                                                   registro, multiple, assays_name)

                write_chain(registro, registro_anterior, cadena, worksheet2)
                write_chain2(registro, registro_anterior, chain3, worksheet2, lengths, variant_chain,
                             assays_name, df_corr)

                registro_anterior = registro

            workbook.close()


def inicio():
    app = QApplication([])  # Crear  aplicacion
    application = MainWindow()  # Llamar a la clase principal
    application.show()  # Mostrar
    # sys.exit(app.exec())  # Salir
    app.exec()


def filtros_final(filename):
    read_file = pd.read_excel(filename, sheet_name='Bad')
    lista_referencias = read_file.Registro
    lista_referencias = lista_referencias.drop_duplicates()

    amplicon = []
    probeok = []
    lista_1 = [0, 1, 2, 3, 4]
    lista_2 = [10, 11, 12, 13, 14, 20, 21, 22, 23, 24, 30, 31, 32, 33, 34, 40, 41, 42, 43, 44]

    chancu_ceros(read_file['Chain'], lista_1, lista_2)
    chancu_ceros(read_file['CHAIN2'], lista_1, lista_2)
    chancu_ceros(read_file['SeqGen'], lista_1, lista_2)

    for i in lista_referencias:
        var_check = read_file.loc[(read_file['Registro'] == i)]

        resta = max(var_check['Inicio']) - min(var_check['Inicio'])

        if len(var_check) == 3:
            if -500 < resta < 500:
                amplicon.append('Ok')
                amplicon.append('Ok')
                amplicon.append('Ok')
            else:
                amplicon.append('No Ok')
                amplicon.append('No Ok')
                amplicon.append('No Ok')

        inicio_forward = var_check.iloc[0]['Inicio']
        inicio_reverse = var_check.iloc[1]['Inicio']
        inicio_probe = var_check.iloc[2]['Inicio']

        if int(inicio_forward) < int(inicio_probe) < int(inicio_reverse) or \
                int(inicio_reverse) < int(inicio_probe) < int(inicio_forward):
            probeok.append('Ok')
            probeok.append('Ok')
            probeok.append('Ok')
        else:
            probeok.append('No Ok')
            probeok.append('No Ok')
            probeok.append('No Ok')

    read_file['Amplicon'] = pd.DataFrame(amplicon)
    read_file['ProbeOk'] = pd.DataFrame(probeok)

    return read_file


def chancu_ceros(row, lista_1, lista_2):
    for w, j in zip(row, range(len(row))):
        if w in lista_1:
            row[j] = '00' + str(w)

        elif w in lista_2:
            row[j] = '0' + str(w)

        else:
            row[j] = str(row[j])
