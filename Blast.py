import os
import time
import pkg_resources


def blast(primers, patogeno, evalue, maxseqs, extrar, extrap):
    """

    :param primers: Contiene los primers que queremos buscar en blast.
    :param patogeno: El nombre del archivo que se va a generar.
    :return:
    """
    tic = time.perf_counter()
    folder = pkg_resources.resource_filename(__name__, "folder")
    # os.chdir('./bin')  # Cambiar al directorio ./bin, donde se va a generar el archivo.
    e_value = " -evalue " + evalue
    init = folder + '\\' + "blastn -task blastn-short "
    db = "-db nt "
    max_seqs = " -max_target_seqs " + maxseqs
    out = "-out " + patogeno
    extra = " -reward " + extrar + " -penalty " + extrap
    form = " -outfmt " + '"6 qseqid sacc sstart send qstart qend mismatch sskingdoms sblastnames sscinames ' \
                         'sstrand evalue qseq sseq" '

    totalquery = init + db + "-query " + primers + e_value + ' -remote ' + out + max_seqs + extra + form
    print(totalquery)

    os.system(totalquery)  # Hacer la busqueda.

    toc = time.perf_counter()
    print(f"Code in {toc - tic:0.4f} seconds")
    print("FINISH BLAST")


