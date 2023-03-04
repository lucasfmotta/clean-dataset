#!/usr/bin/python3

import os
import sys
from Bio import SeqIO
from tqdm import tqdm

def get_files():
    lista_dbs = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.db')]
    lista_fastqs = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.fastq')]
    if not lista_dbs or not lista_fastqs:
        print('No .db or .fastq files found')
        quit()
    make_dbs(lista_dbs)
    return lista_dbs, lista_fastqs

def make_dbs(lista_dbs):
    '''
    genera las db para blast
    '''
    for db in lista_dbs:
        problema = os.system(f'makeblastdb -dbtype nucl -in {db}')
        if problema:
            quit()

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
        problema = os.system(f'blastn -query query.fasta -db {db} -out query.resultado {parametros}')
        if problema:
            quit()
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

def clean_tempfiles():
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

def read_input():
    if '-h' in sys.argv or '-help' in sys.argv:
        show_help()
    if '-U' in sys.argv:
        update()
    if '-query' in sys.argv or '-db' in sys.argv or '-out' in sys.argv:
        print('Wrong input. For help use the -h command')
        quit()
    else:
        parametros_blast = (' ').join(sys.argv[1:])
        return parametros_blast

def show_help():
    print('clean-dataset.py [BLAST parameters]\n  [optional]\nUpdate: -U')
    quit()


parametros_blast = read_input()
lista_dbs, lista_fastqs = get_files()

if not os.path.isdir('Results'):
    os.mkdir('Results')

salida_hit = open('Results/HIT.fastq', 'w')
salida_nohit = open('Results/NO-HIT.fastq', 'w')

n_seqs_total = count_seqs(lista_fastqs)
n_parcial = 0

for fastq in tqdm(lista_fastqs):
    for seq in SeqIO.parse(fastq, 'fastq'):
        n_parcial += 1
        SeqIO.write(seq, 'query.fasta', 'fasta')
        if blast(lista_dbs, parametros_blast): # hubo hit con alguna db -> HIT.fastq
            salida_hit.write(seq.format('fastq'))
        else:
            salida_nohit.write(seq.format('fastq')) # No hubo hit con ninguna db -> NO-HIT.fastq
        porcentaje = n_parcial / n_seqs_total * 100
        print(f'\r{round(porcentaje, 1)}%', end='')

salida_hit.close()
salida_nohit.close()
clean_tempfiles()
print('Goodbye')  # solo para que quede lindo el final