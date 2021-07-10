# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


import random

# import numpy as np
# import scipy.sparse

from faces import Faces
from smeg_matrix import *


file_path = '/Users/ruslanpepa/PycharmProjects/yamabe_flows_tetrahedron/tetrahedron.txt'
VERTEX = 4  # количество вершин в многограннике
EDGES = 6  # количество ребер в многограннике
FACES = 4  # количестов граней в многограннике
TIMES = 1000  # количество шагов по времени
step_time = 0.001 # шаг по времени
list_faces = []  # список, который будет содержать все грани
with open(file_path) as fl_wth_fs:  # выгрузим из файла все номера вершин
    lines = fl_wth_fs.readlines()
for line in lines:  # все номера вершин загоним в списко файлов
    ns_vx = line.rstrip('\n').split('\t')  # получили только числа из каждой строки
    a = int(ns_vx[0])
    b = int(ns_vx[1])
    c = int(ns_vx[2])
    list_faces.append(Faces(a, b, c))
conformal_weights = np.zeros((VERTEX, TIMES), float)  # конформные веса в вершинах
for i in range(0, VERTEX):
    conformal_weights[i, 0] = 1
gauss_curve = adjacency_matrix(list_faces, VERTEX)  # гауссова кривизна в вершинах многогранника
length_matrix = adjacency_matrix(list_faces, VERTEX)  # матрица смежности длин рёбер
trials = 0
while True: # запускаем цикл, который образом создаёт тетраэдр с случайным набором длин рёбер
    # for i in range(0, length_matrix.count_nonzero()):
    #     row, col = length_matrix.nonzero()  # список все индексов в строке, которые
    #     length_matrix[row[i], col[i]] = length_matrix[row[i], col[i]] = random.uniform(1, 10) # здесь присваивается значения различным длинам рёбер
    for i in range(0, VERTEX):
        for j in range(i, VERTEX):
            if length_matrix[i, j ] != 0:
                length_matrix[i, j] = length_matrix[j, i] =  random.uniform(1, 10)
    try:
        trials += 1
        number_of_vertex = len(gauss_curve_calculate(length_matrix)) # проверяем, не схлопнулась ли грань
    except:
        number_of_vertex = 0
    if number_of_vertex == VERTEX: # если все грани в порядке, то завершаем цикл
        print('trials: ',trials)
        break
gauss_curve = gauss_curve_calculate(length_matrix)
# print(np.pi)
print((gauss_curve))
# print(length_matrix.toarray())
# print(сayley_menger_determinant(length_matrix, VERTEX))
for i in range(0, TIMES-1):
    for j in range(0, VERTEX):
        # k1 = k2 = k3 = k4 = .0
        k1 = -(gauss_curve[j] - 2.*np.pi / VERTEX)*conformal_weights[j, i]
        k2 = -(gauss_curve[j] - 2.*np.pi / VERTEX)*(conformal_weights[j, i] + step_time*k1/2.)
        k3 = -(gauss_curve[j] - 2.*np.pi / VERTEX)*(conformal_weights[j, i] + step_time*k2/2.)
        k4 = -(gauss_curve[j] - 2.*np.pi / VERTEX)*(conformal_weights[j, i] + step_time*k3)
        conformal_weights[j, i+1] = conformal_weights[j, i] + (step_time/6.)*(k1 + k2*2. + k3*2. + k4)
    vector_times = conformal_weights[:, i]
    print(vector_times)





# caley_menger = None
# while caley_menger != None:
#     for i in range(0, length_matrix.count_nonzero()):
#         row, col = length_matrix.nonzero()  # список все индексов в строке, которые
#         length_matrix[row[i], col[i]] = random.uniform(1, 10)
#         print('hello, i am here hello, i am here hello, i am here hello, i am here hello, i am here')
# caley_menger = gauss_curve_calculate(length_matrix)
# print(caley_menger)
# print(сayley_menger_determinant(length_matrix, VERTEX))