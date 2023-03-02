#!/usr/bin/python3

import os
import sys

def make_dbs(lista_dbs):
    '''
    genera las db para blast
    '''
    for db in lista_dbs:
        os.system(f'makeblastdb -dbtype nucl -in {db}')

def make_query_fasta(lineas_query):
    '''
    genera fasta con la secuencia query
    '''
    salida = open('query.fasta', 'w')
    salida.write(f'>{lineas_query[1:]}')  # cambia @ por >
    salida.close()

def read_result():
    '''
    lee resultado del blast.
    devuelve False si no hubo hit/True si hubo hit
    '''
    resultado = open('query.resultado')
    for line in resultado:
        if line.find('No hits found') != -1:  # si no hubo hit con esa db
            resultado.close()
            os.remove('query.resultado')
            return False
    resultado.close()
    os.remove('query.resultado')
    return True

def blast(lista_dbs, parametros):
    '''
    corre el blast con la secuencia query y todas las dbs.
    cuando hay hit con una db, corta asi no sigue probando con las otras dbs al pedo
    '''
    for db in lista_dbs:
        os.system(f'blastn -query query.fasta -db {db} -out query.resultado {parametros}')
        if read_result():  # si hubo hit
            os.remove('query.fasta')
            return True
    os.remove('query.fasta')
    return False

def count_seqs(lista_fastqs):
    '''
    devuelve numero total de secuencias
    es solo para calcular el porcentaje completado
    '''
    n = 0
    for fastq in lista_fastqs:
        fastq = open(fastq)
        for line in fastq:
            if line.find('runid=') != -1:
                n += 1
        fastq.close()
    return n

def clean_tempfiles():
    temp_files = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.nhr') or xfile.endswith('nin') or xfile.endswith('nsq')]
    for temp_file in temp_files:
        os.remove(temp_file)

def show_help():
    print('clean-dataset.py [BLAST parameters]\n  [optional]')
    quit()

parametros_blast = (' ').join(sys.argv[1:])
if parametros_blast == '-h' or parametros_blast == '-help':
    show_help()

lista_dbs = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.db')]
make_dbs(lista_dbs)

lista_fastqs = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.fastq')]

ids_secuencias_basura = [] # aca se juntan los ids de todas las secuencias que dieron hit con alguna db
n_seqs_total = count_seqs(lista_fastqs)
n_parcial = 0
for fastq in lista_fastqs:
    fastq = open(fastq)
    read_seq = False
    for line in fastq:
        if line.find('runid=') != -1:
            n_parcial += 1
            lineas_query = line
            id_query = line[line.find('read=')+5:line.find(' ', line.find('read='))] # Busca el identificador de la secuencia basado en el número de read
            read_seq = True
            continue
        if read_seq:
            lineas_query += line
            read_seq = False
            make_query_fasta(lineas_query)
            if blast(lista_dbs, parametros_blast): # hubo hit con alguna db, la secuencia es basura
                ids_secuencias_basura.append(id_query)
            porcentaje = n_parcial / n_seqs_total * 100
            print(f'\r{round(porcentaje, 2)}%', end='')
    fastq.close()

dataset_limpio = open('DATASET.fastq', 'w')

for fastq in lista_fastqs:
    fastq = open(fastq)
    for line in fastq:
        if line.find('runid=') != -1:
            id_query = line[line.find('read=') + 5 : line.find(' ', line.find('read='))] # Busca el identificador de la secuencia basado en el número de read
            if id_query in ids_secuencias_basura:
                write = False
                ids_secuencias_basura.remove(id_query) # saca el id de la lista asi las proximas busquedas son mas rapidas
                continue
            else:
                write = True
        if write:
            dataset_limpio.write(line)
    fastq.close()

dataset_limpio.close()
clean_tempfiles()
print('Goodbye')  # solo para que quede lindo el final