# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


import random
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
# import numpy as np
# import pandas as pd
# import numpy as np
# import scipy.sparse

from faces import Faces
from smeg_matrix import *

file_path = '/Users/ruslanpepa/PycharmProjects/yamabe_flows_tetrahedron/tetrahedron.txt'
VERTEX = 4  # количество вершин в многограннике
EDGES = 6  # количество ребер в многограннике
FACES = 4  # количестов граней в многограннике
TIMES = 1000  # количество шагов по времени
step_time = 0.001  # шаг по времени
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
gauss_curvature = np.zeros((VERTEX, TIMES), float) #
for i in range(0, VERTEX):
    conformal_weights[i, 0] = 1
gauss_curve = adjacency_matrix(list_faces, VERTEX)  # гауссова кривизна в вершинах многогранника
length_matrix = adjacency_matrix(list_faces, VERTEX)  # матрица смежности длин рёбер
trials = 0
while True:  # запускаем цикл, который образом создаёт тетраэдр с случайным набором длин рёбер
    # for i in range(0, length_matrix.count_nonzero()):
    #     row, col = length_matrix.nonzero()  # список все индексов в строке, которые
    #     length_matrix[row[i], col[i]] = length_matrix[row[i], col[i]] = random.uniform(1, 10)
    for i in range(0, VERTEX):
        for j in range(i, VERTEX):
            if length_matrix[i, j] != 0:
                length_matrix[i, j] = length_matrix[j, i] = 2.#random.uniform(1, 10)
    number_of_vertex = 0
    try:
        trials += 1
        number_of_vertex = len(gauss_curve_calculate(length_matrix))  # проверяем, не схлопнулась ли грань
    except Exception:
        continue
    if number_of_vertex == VERTEX and keyle_menger_det(length_matrix, VERTEX) > 0:  # завершаем цикл
        print('trials: ', trials)
        break
gauss_curve = gauss_curve_calculate(length_matrix)
print('gauss_curve: hello', gauss_curve)
for i in range(0, VERTEX):
    gauss_curvature[0,i] = gauss_curve[i]
    print(gauss_curvature[0, i])
# print(type(gauss_curve[1]))
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
    length_matrix = get_length(length_matrix, conformal_weights[:, i+1])  # Пересчитываем все длины сторон
    gauss_curve = gauss_curve_calculate(length_matrix)  # Пересчитываем все значения кривизн в вершинах тетраэдра
    if keyle_menger_det(length_matrix, VERTEX) < 0:
        break
    for j in range(0, VERTEX):
        gauss_curvature[j, i] = gauss_curve[j]
    print('gauss curve:', gauss_curve)
    # print(length_matrix.todense())
    # print('number of iteration: ', i)

# caley_menger = None
# while caley_menger != None:
#     for i in range(0, length_matrix.count_nonzero()):
#         row, col = length_matrix.nonzero()  # список все индексов в строке, которые
#         length_matrix[row[i], col[i]] = random.uniform(1, 10)
#         print('hello, i am here hello, i am here hello, i am here hello, i am here hello, i am here')
# caley_menger = gauss_curve_calculate(length_matrix)
# print(caley_menger)
# print(conformal_weights[1,:])
# plt.plot(gauss_curvature[0, 0:-2])
# plt.plot(gauss_curvature[1, 0:-2])
# plt.plot(gauss_curvature[2, 0:-2])
plt.plot(gauss_curvature[3, 0:-2])
# plt.plot([25, 20, 15, 12, 9, 4, 2])
plt.savefig('gauss_curve№3.png')
plt.show()
# import numpy as np
# from matplotlib import pyplot as plt


plt.style.use('seaborn-pastel')

fig = plt.figure()
ax = plt.axes(xlim=(0, step_time*TIMES), ylim=(0, 2*np.pi))
line, = ax.plot([], [], lw=3)


def init():
    line.set_data([], [])
    return line,


def animate(i):
    # x = np.linspace(0, 4, 1000)
    # y = np.sin(2 * np.pi * (x - 0.01 * i))
    x =np.linspace(0, step_time*TIMES, TIMES)
    y = gauss_curvature[0, i]
    line.set_data(x, y)
    return line,


anim = FuncAnimation(fig, animate, init_func=init,
                     frames=200, interval=20, blit=True)
# plt.show()
anim.save('sine_wave.gif', writer='pillow')
