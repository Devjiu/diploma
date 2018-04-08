# -*- coding: utf-8 -*-
import csv
import time


# Эта функция описывает процедуру записи результатов расчетов в файл
def prp_atm(atoms, k):
    with open('atoms_out' + str(k) + '.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        header = ['ID', 'x', 'y', 'z', 'sigma', 'epsilon', 'charge', 'Rvdw', 'volume', 'AA', 'PDB']
        csvwriter.writerow(header)
        for i in atoms:
            csvwriter.writerow(i)


# Эта функция нужна, чтобы следить за временем работы программы
def millis():
    return int(round(time.time() * 1000))
