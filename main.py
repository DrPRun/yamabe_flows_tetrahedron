# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


import random
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from faces import Faces
from smeg_matrix import *

file_path = '/Users/ruslanpepa/PycharmProjects/yamabe_flows_tetrahedron/tetrahedron.txt'
VERTEX = 4  # количество вершин в многограннике
EDGES = 6  # количество ребер в многограннике
FACES = 4  # количестов граней в многограннике
TIMES = 500 # количество шагов по времени
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
conformal_weights = np.ones((VERTEX, TIMES), float)  # конформные веса в вершинах
gauss_curvature = np.zeros((VERTEX, TIMES), float) # гауссова кривизна в начальный момент времени
length_of_tetrahedron = np.zeros((EDGES, TIMES), float) # экспериментальная матрица для отображения длин рёбер
kayli_manger = np.zeros((FACES, TIMES), float) # массив, который будет содержать значения определителей Кэлли-Менгера на грани
# gauss_curve = adjacency_matrix(list_faces, VERTEX)  # гауссова кривизна в вершинах многогранника
#print('gauss_curve:', gauss_curve)
# print(gauss_curve.todense())
length_matrix = adjacency_matrix(list_faces, VERTEX)  # матрица смежности длин рёбер
print(length_matrix.todense())
trials = 0
while True:  # запускаем цикл, который образом создаёт тетраэдр с случайным набором длин рёбер
    # for i in range(0, length_matrix.count_nonzero()):
    #     row, col = length_matrix.nonzero()  # список все индексов в строке, которые
    #     length_matrix[row[i], col[i]] = length_matrix[row[i], col[i]] = random.uniform(1, 10)
    for i in range(0, VERTEX):
        for j in range(i, VERTEX):
            if length_matrix[i, j] != 0:
                length_matrix[i, j] = length_matrix[j, i] = random.uniform(1, 10)
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
    gauss_curvature[0, i] = gauss_curve[i]
    print(gauss_curvature[0, i])
# нижние три строки нужны для того, чтобы попробовать анимировать эволюцию грани
length_of_tetrahedron[0, 0] = length_matrix[0, 1]
length_of_tetrahedron[0, 1] = length_matrix[1, 2]
length_of_tetrahedron[0, 2] = length_matrix[0, 2]

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
        row = [0, 1, 2]
        col = [0, 1, 2]
        data = [a, b, c]
        space_r = np.array(row)
        space_c = np.array(col)
        space_d = np.array(data)
        fases_len_matrix = sparse.coo_matrix((space_d, (space_r, space_c)), shape=(3, 3)).tocsc()
        half_perimetr = (a+b+c)/2.
        kl_mng = np.sqrt(half_perimetr*(half_perimetr - a)*(half_perimetr - b)*(half_perimetr - c))
        kayli_manger[times, i+1] = kl_mng
        times += 1
        print('keyli_menger:',kl_mng)



    length_matrix = get_length(length_matrix, conformal_weights[:, i+1])  # Пересчитываем все длины сторон
    gauss_curve = gauss_curve_calculate(length_matrix)  # Пересчитываем все значения кривизн в вершинах тетраэдра
    if keyle_menger_det(length_matrix, VERTEX) <= 0:
        break
    for j in range(0, VERTEX):
        gauss_curvature[j, i+1] = gauss_curve[j]
    print('gauss curve:', gauss_curve)
    length_of_tetrahedron[0, i+1] = length_matrix[0, 1]
    length_of_tetrahedron[1, i+1] = length_matrix[1, 2]
    length_of_tetrahedron[2, i+1] = length_matrix[0, 2]

# print(kayli_manger)
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
# plt.plot(gauss_curvature[3, 0:-2])
# massiv_sum = sum(gauss_curvature[:, 0:-2])
plt.plot(kayli_manger[0, 0:-2])
plt.plot(kayli_manger[1, 0:-2])
plt.plot(kayli_manger[2, 0:-2])
plt.plot(kayli_manger[3, 0:-2])
# plt.plot(massiv_sum)
# plt.plot(length_matrix[1,2])
# plt.plot([25, 20, 15, 12, 9, 4, 2])
# plt.savefig('gauss_curve№3.png')
plt.show()

# a = length_of_tetrahedron[0, 1]
# b = length_of_tetrahedron[1, 1]
# c = length_of_tetrahedron[2, 1]
# cosine = (a**2 + c**2 - b**2)/(2*a*c)
# sines = np.sqrt(1 - cosine**2)
# fig = plt.figure(figsize=(10, 10))
# ax = plt.axes([0, 0, 1, 1])
# triangle1 = mpatches.Polygon(np.array([[0,0],[a*cosine,a*sines],[c,0]]), fc="blue")
# ax.add_artist(triangle1)
# ax.set_xlim(-3, 3)
# ax.set_ylim(-3, 3)
# plt.show()


# fig = plt.figure()
# ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
# line, = ax.plot([], [], lw=2)
#
# # initialization function: plot the background of each frame
# def init():
#     line.set_data([], [])
#     return line,
#
# # animation function.  This is called sequentially
# def animate(i):
#     a = length_of_tetrahedron[0, i+1]
#     b = length_of_tetrahedron[1, i+1]
#     c = length_of_tetrahedron[2, i+1]
#     cosine = (a**2 + c**2 - b**2)/(2*a*c)
#     sines = np.sqrt(1 - cosine**2)
#     x = [0.0, a*cosine,  c]
#     y = [0.0, a*sines, 0]
#     line.set_data(x, y)
#     return line,
#
# # call the animator.  blit=True means only re-draw the parts that have changed.
# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                frames=200, interval=20, blit=True)
#
# # save the animation as an mp4.  This requires ffmpeg or mencoder to be
# # installed.  The extra_args ensure that the x264 codec is used, so that
# # the video can be embedded in html5.  You may need to adjust this for
# # your system: for more information, see
# # http://matplotlib.sourceforge.net/api/animation_api.html
# anim.save('sine_wave.gif', writer='pillow')
#
# plt.show()
