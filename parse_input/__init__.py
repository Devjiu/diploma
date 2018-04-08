# -*- coding: utf-8 -*-
import os, sys

nb_dir = os.path.split(os.getcwd())[0]
if nb_dir not in sys.path:
    sys.path.append(nb_dir)


# Функция переводит считываемый .txt файл (аргумент str)
# без первой строчки (не содержащей данных) в список python
# (обект типа list); аргумент m отвечает за то,
# сколько в начале идет переменных int-типа; аргумент M был введен,
# чтобы не считывать не нужную информацию
def rd(name_of_file, m=100, M=100):
    def path_to_resources(filename):
        script_dir = os.path.split(os.getcwd())
        # print(script_dir)
        abs_file_path = os.path.abspath(os.path.join(script_dir[0], "data", name_of_file))
        return abs_file_path
    absolute_path = path_to_resources(name_of_file)
    with open(absolute_path) as f:
        poly_shape = []
        k = 0
        for line in f:
            if k > 0:
                line = line.strip()
                line = line.split(",")
                newline = line[:min(M, len(line))]  # to deal with blank
                if newline:  # lines (ie skip them)
                    newline = [float(i) for i in newline]
                    l = len(newline)
                    ml = min(l, m)
                    for i in range(0, ml):
                        newline[i] = int(newline[i])
                    poly_shape.append(newline + line[min(M, len(newline)):])
            k = k + 1
    return poly_shape