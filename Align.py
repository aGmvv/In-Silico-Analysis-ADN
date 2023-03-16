import os
import csv
import sys
import time
import ali_ui
import operator
import xlsxwriter
import pandas as pd
from Bio import SeqIO
from itertools import zip_longest
from reverso_complementario import main
from PyQt5.QtWidgets import QApplication, QComboBox, QMainWindow, QTabWidget, QFileDialog


def classifier2(row):
    if (int(row.Fstart) < int(row.Pstart) < int(row.Rstart)) or (int(row.Rstart) < int(row.Pstart) < int(row.Fstart)):
        row['ProbeOk'] = 'Yes'
    else:
        row['ProbeOk'] = 'No'
    return row


def process(datamerged):
    print("RTETETEa")
    interF = datamerged.query('Oligo == "FORWARD"')
    dmF = interF.groupby(['Referencia']).Oligo.count()
    dmF = pd.DataFrame(dmF)
    dmF = dmF.reset_index()
    interR = datamerged.query('Oligo == "REVERSE"')
    dmR = interR.groupby(['Referencia']).Oligo.count()
    dmR = pd.DataFrame(dmR)
    dmR = dmR.reset_index()
    interP = datamerged.query('Oligo == "PROBE"')
    dmP = interP.groupby(['Referencia']).Oligo.count()
    dmP = pd.DataFrame(dmP)
    dmP = dmP.reset_index()
    regs = datamerged.Referencia.unique()
    regs = pd.DataFrame(regs, columns=['Referencia'])
    newdata = regs.merge(dmF, on='Referencia', how='left')
    newdata = newdata.merge(dmR, on='Referencia', how='left')
    newdata = newdata.merge(dmP, on='Referencia', how='left')
    newdata.columns = ['Referencia', 'FORWARD', 'REVERSE', 'PROBE']
    newdata2 = newdata[(newdata.FORWARD >= 1) & (newdata.REVERSE >= 1) & (newdata.PROBE >= 1)]
    ls = list(newdata2.Referencia)
    finalls = pd.DataFrame(ls, columns=['Referencia'])
    newdata3 = datamerged.query('Referencia in @ls')
    nF = newdata3[newdata3['Oligo'] == 'FORWARD']
    nR = newdata3[newdata3['Oligo'] == 'REVERSE']
    nP = newdata3[newdata3['Oligo'] == 'PROBE']
    nF = nF[['Referencia', 'Rstart']]
    nR = nR[['Referencia', 'Rstart']]
    nP = nP[['Referencia', 'Rstart']]
    nF = nF.drop_duplicates(subset=['Referencia'])
    nR = nR.drop_duplicates(subset=['Referencia'])
    nP = nP.drop_duplicates(subset=['Referencia'])
    return newdata3, finalls, nF, nR, nP


def send_mail(to, file):
    pass


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
        # worksheet2.write(x, 18, reason)
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
            # print('MAX1= ', max(first3), ' MAX2= ', max(second3))
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
    mmax = 1
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
        if tg[1] == '':
            if oliname == '':
                if seq == '':
                    variant = ''
        if tg[1] == '' and tg[5] == '':
            if oliname == '':
                if seq == '':
                    variant = ''
                if seq == '':
                    variant = ''
                if seq == '':
                    variant = ''
            if oliname == '':
                if seq == '':
                    variant = ''
                if seq == '':
                    variant = ''
        if tg[1] == '' and tg[5] == '':
            if oliname == '':
                if seq == '':
                    variant = ''
        if tg[5] == '':
            if oliname == '':
                if seq == '':
                    variant = ''
        if tg[5] == '':
            if oliname == '':
                if seq == '':
                    variant = ''
                if seq == '':
                    variant = ''
                if seq == '':
                    variant = ''
        if tg[5] == '':
            if oliname == '':
                if seq == '':
                    variant = ''
        if tg[1] == '' and tg[5] == '':
            if oliname == '':
                if seq == '':
                    variant = ''
            if oliname == '':
                if seq == '':
                    variant = ''
    return variant


def create_file(fw, rv, sd):
    file = open("./set2.fsa", "w")
    file.write(">FORWARD\n")
    file.write(fw + "\n")
    file.write(">REVERSE\n")
    file.write(rv + "\n")
    file.write(">PROBE\n")
    file.write(sd + "\n")
    file.close()
    print("Archivo Creado")


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


def get_threads():
    core = int(os.environ['NUMBER_OF_PROCESSORS'])
    if core > 1:
        core = core - 1
    print("Numero de cores: " + str(core))
    return core


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
        self.ui.pushButton_2.clicked.connect(self.resume)
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

    def xmlparser(self):
        from Bio.Blast import NCBIXML
        print("PARSEANDO" + self.ui.lineFile_2.text())
        result_handle = open(self.ui.lineFile_2.text())
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            print(blast_record.query)
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    print("****Alignment****")
                    print("SeqName", )
                    print("sequence:", alignment.title)
                    print("length:", alignment.length)
                    print("e value:", hsp.expect)
                    print(hsp.query[0:75] + "...")
                    print(hsp.match[0:75] + "...")
                    print(hsp.sbjct[0:75] + "...")

    def blast(self):
        tic = time.perf_counter()
        fw = self.ui.lineEditBlFW.text()
        rv = self.ui.lineEditBlRV.text()
        sd = self.ui.lineEditBlSD.text()
        name = self.ui.lineNombreOut.text()
        e_value = " -evalue " + str(self.ui.lineMaxSeqs_2.text())
        remote = ""
        cores = " -num_threads " + str(get_threads())
        if (self.ui.checkBoxRemote.checkState()) == 2:
            remote = " -remote "
            cores = ""
        init = "blastn -task blastn-short"
        db = " -db " + self.ui.lineDb.text()
        max_seqs = " -max_target_seqs " + self.ui.lineMaxSeqs.text()
        out = " -out " + name
        print("Numero de cores: " + str(cores))
        if self.ui.comboBox.currentIndex() == 6:
            form = " -outfmt " + '"6 qseqid sacc sstart send qstart qend mismatch sskingdoms sblastnames sscinames ' \
                                 'sstrand evalue" '
        else:
            form = " -outfmt " + str(self.ui.comboBox.currentIndex())
        create_file(fw, rv, sd)
        print(form)
        print(type(form))

        totalquery = init + db + " -query set2.fsa " + e_value + cores + remote + out + max_seqs + form
        print(totalquery)
        os.system('cmd /c' + totalquery)

        # Emepezamos el resumen de pandas
        if self.ui.comboBox.currentIndex() == 6:
            df = pd.read_csv(name, sep="\t", header=None)
            records = list(SeqIO.parse("set2.fsa", "fasta"))
            olidf = pd.DataFrame({'LENGTH': [len(records[0].seq), len(records[1].seq), len(records[2].seq)],
                                  'Oligo': ['FORWARD', 'REVERSE', 'PROBE']})
            df.columns = ['Oligo', 'Referencia', 'Rstart', 'Rend', 'Ostart', 'Oend', 'Mismatch', 'Kingdom', 'Group',
                          'Group2', 'Strand', 'Evalue']

            df['maxL'] = df.Oend - df.Ostart + 1
            datamerged = df.merge(olidf, on='Oligo')
            datamerged['MM'] = datamerged.LENGTH - datamerged.maxL
            newdata3, finalls, nF, nR, nP = process(datamerged)

            final = finalls.merge(nF[['Referencia', 'Rstart']], on='Referencia')
            final = final.merge(nP[['Referencia', 'Rstart']], on='Referencia')
            final = final.merge(nR[['Referencia', 'Rstart']], on='Referencia')
            final.columns = ['Referencia', 'Fstart', 'Pstart', 'Rstart']
            final = final.merge(newdata3[['Referencia', 'Group', 'Kingdom', 'Group2']], on='Referencia')
            final = final.drop_duplicates()

            final['resta'] = final.Rstart - final.Fstart
            final = final.fillna('Indeterminado')
            final2 = final.groupby(['Kingdom', 'Group', 'Group2', ]).Referencia.count()
            final = final.apply(classifier2, axis='columns')
            print(final)
            final3 = final[(final.ProbeOk == 'Yes')].groupby(['Kingdom', 'Group', 'Group2', ]).Referencia.count()
            print("RTETETEa3")
            ex_name = name + "_cross.xlsx"
            print(datamerged)
            datamerged.to_excel(ex_name, sheet_name='raw')
            with pd.ExcelWriter(ex_name, mode='a', engine="openpyxl") as writer:
                final.to_excel(writer, sheet_name='data')
                final2.to_excel(writer, sheet_name='resume')
                final3.to_excel(writer, sheet_name='ProbeOk')
        toc = time.perf_counter()
        print(f"Code in {toc - tic:0.4f} seconds")
        print("FINISH BLAST")

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

        if s == "":
            data = []

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
        if hasattr(sys, '_MEIPASS'):
            path = sys._MEIPASS
        else:
            path = os.getcwd()

        where = os.path.join(path, "Sources/correla2.xlsx")
        df_corr = pd.read_excel(where, dtype=str, engine='openpyxl')

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
        # send_mail(self.ui.lineEmail.text(), self.ui.lineFile.text())
        toc = time.perf_counter()
        print(f"Code in {toc - tic:0.4f} seconds")

    def ali(self, ford, reve, sond, ford2, reve2, sond2, ford3, reve3, sond3, ford4, reve4, sond4, file2, target1,
            target2, target3, target4, ndiana, df_corr):
        from Bio import SeqIO
        from Bio import Align

        one_and_wrong = []

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

                    # --------------------------------------------------------------------------------------------- #
                    #                              Alinea y asigna numero de gravedad                               #
                    # --------------------------------------------------------------------------------------------- #
                    OlAlign = aligner.align(x, records[i].seq.upper())
                    OlAlign_reverso = aligner.align(main(x), records[i].seq.upper())

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
