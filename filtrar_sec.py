def filtrar_archivo(filename):
    print('Creando archivo...')

    with open(filename, 'r') as f_in, \
            open('archivo_filtrado.txt', 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                f_out.write(line)

    print('Archivo creado')
