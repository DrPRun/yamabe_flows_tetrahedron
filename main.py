# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


import random

import matplotlib.pyplot as plt

from faces import Faces
from smeg_matrix import *
from animation import *

file_path = '/Users/ruslanpepa/PycharmProjects/yamabe_flows_tetrahedron/tetrahedron.txt'
VERTEX = 4  # количество вершин в многограннике
EDGES = 6  # количество ребер в многограннике
FACES = 4  # количестов граней в многограннике
TIMES = 100 # количество шагов по времени
step_time = 0.01  # шаг по времени
list_faces = []  # список, который будет содержать все грани
with open(file_path) as fl_wth_fs:  # выгрузим из файла все номера вершин
    lines = fl_wth_fs.readlines()
for line in lines:  # все номера вершин загоним в списко файлов
    ns_vx = line.rstrip('\n').split('\t')  # получили только числа из каждой строки
    a = int(ns_vx[0])
    b = int(ns_vx[1])
    c = int(ns_vx[2])
    list_faces.append(Faces(a, b, c))
conformal_weights = np.ones((VERTEX, TIMES), float)  # конформные веса в вершинах
gauss_curvature = np.zeros((VERTEX, TIMES), float) # гауссова кривизна в начальный момент времени
length_of_tetrahedron = np.ones((EDGES, TIMES), float) # экспериментальная матрица для отображения длин рёбер
kayli_manger = np.zeros((FACES, TIMES), float) # массив, который будет содержать значения определителей Кэлли-Менгера на грани
length_matrix = adjacency_matrix(list_faces, VERTEX)  # матрица смежности длин рёбер
# print('len(length_of_tetrahedron[1, :]):', len(length_of_tetrahedron[1, :]))
# print(length_matrix.todense())
trials = 0
while True:  # запускаем цикл, который образом создаёт тетраэдр с случайным набором длин рёбер
    for i in range(0, VERTEX):
        for j in range(i, VERTEX):
            if length_matrix[i, j] != 0:
                length_matrix[i, j] = length_matrix[j, i] = random.uniform(8, 9)
    trials += 1
    print(trials)
    if len(gauss_curve_calculate(length_matrix)) != VERTEX or keyle_menger_det(length_matrix, VERTEX) <= 0:
        continue
    else:
        break
gauss_curve = gauss_curve_calculate(length_matrix)
# print('gauss_curve: hello', gauss_curve)
for i in range(0, VERTEX):
    gauss_curvature[0, i] = gauss_curve[i]
    print(gauss_curvature[0, i])
# нижние три строки нужны для того, чтобы попробовать анимировать эволюцию грани
length_of_tetrahedron[0, 0] = length_matrix[0, 1]
length_of_tetrahedron[1, 0] = length_matrix[1, 2]
length_of_tetrahedron[2, 0] = length_matrix[0, 2]

# print(type(lgauss_curve[1]))
# df = pd.DataFrame({'gauss_VTX_№0': gauss_curve[0], 'gauss_VTX_№1': gauss_curve[1], 'gauss_VTX_№2': gauss_curve[2], 'gauss_VTX_№3': gauss_curve[3]})
# print(df)

for i in range(0, TIMES - 1):
    for j in range(0, VERTEX):
        # k1 = k2 = k3 = k4 = .0
        k0 = -(gauss_curve[j] - 4. * np.pi / VERTEX) * conformal_weights[j, i]
        k1 = -(gauss_curve[j] - 4. * np.pi / VERTEX) * (conformal_weights[j, i] + step_time * k0 / 2.)
        k2 = -(gauss_curve[j] - 4. * np.pi / VERTEX) * (conformal_weights[j, i] + step_time * k1 / 2.)
        k3 = -(gauss_curve[j] - 4. * np.pi / VERTEX) * (conformal_weights[j, i] + step_time * k2)
        # print('(step_time / 6.) * (k0 + k1 * 2. + k2 * 2. + k3):', (step_time / 6.) * (k0 + k1 * 2. + k2 * 2. + k3))
        conformal_weights[j, i + 1] = conformal_weights[j, i] + (step_time / 6.) * (k0 + k1 * 2. + k2 * 2. + k3)
    # vector_times = conformal_weights[:, i]
    # print('kaly_menger:', keyle_menger_det(length_matrix, VERTEX))
    times = 0
    for fs in list_faces:
        a = length_matrix[fs[0], fs[1]]
        b = length_matrix[fs[1], fs[2]]
        c = length_matrix[fs[2], fs[0]]
        row = [0, 0, 1, 1, 2, 2]
        col = [1, 2, 0, 2, 0, 1]
        data = [a, c, a, b, c, b]
        space_r = np.array(row)
        space_c = np.array(col)
        space_d = np.array(data)
        fases_len_matrix = sparse.coo_matrix((space_d, (space_r, space_c)), shape=(3, 3)).tocsc()
        # print(fases_len_matrix.todense())
        half_perimetr = (a+b+c)/2.
        if half_perimetr*(half_perimetr - a)*(half_perimetr - b)*(half_perimetr - c) <= 0:
            break
        else:
            kl_mng = np.sqrt(half_perimetr*(half_perimetr - a)*(half_perimetr - b)*(half_perimetr - c))
        kayli_manger[times, i+1] = kl_mng
        times += 1
        # print('keyli_menger:',kl_mng)



    length_matrix = get_length(length_matrix, conformal_weights[:, i+1])  # Пересчитываем все длины сторон
    gauss_curve = gauss_curve_calculate(length_matrix)  # Пересчитываем все значения кривизн в вершинах тетраэдра
    if len(gauss_curve) != VERTEX:
        break
    if keyle_menger_det(length_matrix, VERTEX) <= 0 :
        break
    for j in range(0, VERTEX):
        gauss_curvature[j, i+1] = gauss_curve[j]
    # print('gauss curve:', gauss_curve)
    length_of_tetrahedron[0,i+1 ] = length_matrix[0, 1]
    length_of_tetrahedron[1, i+1] = length_matrix[1, 2]
    length_of_tetrahedron[2, i+1] = length_matrix[0, 2]


plt.plot(kayli_manger[0, 0:-2])
plt.plot(kayli_manger[1, 0:-2])
plt.plot(kayli_manger[2, 0:-2])
# plt.plot(massiv_sum)
# plt.plot(length_matrix[1,2])
# plt.plot(sum(gauss_curvature[:, 0:-2]))
plt.show()