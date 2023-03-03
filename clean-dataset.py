#!/usr/bin/python3

import os
import sys
from Bio import SeqIO

def make_dbs(lista_dbs):
    '''
    genera las db para blast
    '''
    for db in lista_dbs:
        os.system(f'makeblastdb -dbtype nucl -in {db}')

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
            n += 1
        fastq.close()
    return n/4

def clean_tempfiles(multiple_tempfiles):
    if multiple_tempfiles:
        temp_files = ['temp_files/'+xfile for xfile in os.listdir(os.curdir+'/temp_files') if xfile.endswith('.fastq')]
        for temp_file in temp_files:
            os.remove(temp_file)
        os.rmdir('temp_files')
    temp_files = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.nhr') or xfile.endswith('.nin') or xfile.endswith('.nsq')]
    for temp_file in temp_files:
        os.remove(temp_file)

def update():
    '''
    actualiza el script a su ultima version
    '''
    os.system('sudo wget https://raw.githubusercontent.com/lucasfmotta/clean-dataset/main/clean-dataset.py -O /usr/local/bin/clean-dataset')
    os.system('sudo chmod a+rx /usr/local/bin/clean-dataset')
    quit()

def show_help():
    print('clean-dataset.py [BLAST parameters]\n  [optional]\nUpdate: -U')
    quit()

parametros_blast = (' ').join(sys.argv[1:])
if parametros_blast == '-h' or parametros_blast == '-help':
    show_help()
if parametros_blast == '-U':
    update()

lista_dbs = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.db')]
make_dbs(lista_dbs)

lista_fastqs = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.fastq')]

seqs_nodescartadas = []
multiple_tempfiles = False
n_seqs_total = count_seqs(lista_fastqs)
n_parcial = 0

for fastq in lista_fastqs:
    for seq in SeqIO.parse(fastq, 'fastq'):
        n_parcial += 1
        SeqIO.write(seq, 'query.fasta', 'fasta')
        if not blast(lista_dbs, parametros_blast): # no hubo hit con alguna db, la secuencia no se descarta
            seqs_nodescartadas.append(seq)
        porcentaje = n_parcial / n_seqs_total * 100
        print(f'\r{round(porcentaje, 1)}%', end='')

        if len(seqs_nodescartadas) == 4000:
            if not multiple_tempfiles: #solo con el primer tempfile
                os.mkdir('temp_files')
            SeqIO.write(seqs_nodescartadas, f'temp_files/temp{n_parcial}.fastq', 'fastq')
            seqs_nodescartadas = []
            multiple_tempfiles = True

os.mkdir('Results')
if multiple_tempfiles:
    if len(seqs_nodescartadas) > 0: #por si quedó alguna secuencia colgada
        SeqIO.write(seqs_nodescartadas, f'temp_files/temp{n_parcial}.fastq', 'fastq')
    salida = open('Results/DATASET-LIMPIO.fastq', 'w')
    temp_files = ['temp_files/'+xfile for xfile in os.listdir(os.curdir+'/temp_files') if xfile.endswith('.fastq')]
    for temp_file in temp_files:
        temp_file = open(temp_file)
        for line in temp_file:
            salida.write(line)
    salida.close()
else:
    SeqIO.write(seqs_nodescartadas, 'Results/DATASET-LIMPIO.fastq', 'fastq')

clean_tempfiles(multiple_tempfiles)
print('\nGoodbye')  # solo para que quede lindo el final