from pandas import read_csv
import numpy as np


def area_by_shoelace(coords):
    x, y = zip(*coords)
    return 0.5 * np.abs(np.dot(x[:-1], y[1:]) + x[-1]*y[0] -
                        np.dot(y[:-1], x[1:]) - y[-1]*x[0])


def dist(point1, point2):
    return ((point2[0]-point1[0])**2+(point2[1]-point1[1])**2)**0.5


def find_splitting_point(triangle, area):
    p1, p2, p3 = triangle
    k = p2[0]*(p1[1]-p3[1]) + p1[0]*(p3[1]-p2[1]) + p3[0]*(p2[1]-p1[1])
    x = (2*area*(p3[0]-p2[0])+p2[0]*k) / k
    y = (p2[1]-p3[1])*(x-p3[0])/(p2[0]-p3[0]) + p3[1]
    return [x, y]


def read_polygon_data(filename):
    file = read_csv(filename, sep=';', decimal=',')
    return file.astype({"X": float, "Y": float})


def get_polygon_coords(dataframe, inverse=False):
    columns = ["X", "Y"]
    if inverse:
        columns.reverse()
    coords = dataframe[columns].values.tolist()
    return coords


def split_polygon(coords):
    total_polygon_area = area_by_shoelace(coords)
    half_area = total_polygon_area*0.5
    triangles = [[coords[0], *coords[i:i+2]] for i in range(len(coords)-2)]
    current_area = 0

    triangle_areas = {i: area_by_shoelace(
        triangle) for i, triangle in enumerate(triangles)}

    for index, triangle in enumerate(triangles):
        current_area += triangle_areas.get(index)
        if current_area >= half_area:
            break

    area = triangle_areas.get(index) - (current_area-half_area)
    splitting_point = find_splitting_point(
        triangle, area)
    split_line = [coords[0], splitting_point]
    coords.insert(index+1, splitting_point)
    return [coords, split_line]


def square_xy(file):
    z = float(0)
    mn1 = float(0)
    mn2 = float(0)
    for x in range(len(file)-1):
        if x == 0:
            xi = float(file.iloc[x, 1])
            yi_1 = float(file.iloc[x+1, 2])
            yi_i1 = float(file.iloc[len(file)-2, 2])
        elif x == len(file) - 2:
            xi = float(file.iloc[x, 1])
            yi_1 = float(file.iloc[0, 2])
            yi_i1 = float(file.iloc[x-1, 2])
        else:
            xi = float(file.iloc[x, 1])
            yi_1 = float(file.iloc[x+1, 2])
            yi_i1 = float(file.iloc[x-1, 2])
        mn1 += xi*yi_1
        mn2 += xi*yi_i1
        if x == len(file) - 2:
            z = float((mn1 - mn2)/2)
            return z


def transform_coordinates(coords, start, stop):
    basis = dist(start, stop)
    cos, sin = np.divide(np.subtract(stop, start), basis)
    matrix = [[cos, sin], [-sin, cos]]
    transformed = map(lambda point: np.matmul(
        matrix, np.subtract(point, start)).tolist(), coords)
    return list(transformed)


# def transform_coordinates(coords, start, stop):
#     x0, y0 = start
#     x1, y1 = stop
#     basis = dist(start, stop)
#     sinA = (y1 - y0)/basis
#     cosA = (x1 - x0)/basis
#     return [[(x - x0)*cosA + (y - y0)*sinA, -
#              (x - x0)*sinA + (y - y0)*cosA] for x, y in coords]


def smooth_polygon(file, status, output_dir):
    start_row = file[file.STATUS.eq("start")]
    stop_row = file[file.STATUS.eq("stop")]
    start = start_row[["X", "Y"]].values[0]
    stop = stop_row[["X", "Y"]].values[0]
    xn, yn = start
    xn_j, yn_j = stop
    # data for new square
    n_s = file[file.STATUS.ne('+')]

    # main square
    d_s = square_xy(file) - square_xy(n_s)  # difference of squares

    # def - для локальной СК, чтобы сразу считать
    basis = dist([xn_j, yn_j], [xn, yn])
    cosA = (yn_j - yn)/basis
    sinA = (xn_j - xn)/basis
    H = 2*d_s/basis

    # Точки сглаживания в локальной СК.
    # По значению Y_l_NEW можно понять, где входит H.
    part = file.dropna()
    len_part = len(part)
    dict_x = []
    dict_y = []
    for l in range(len_part):
        x_l = float(part.iloc[l, 1])
        y_l = float(part.iloc[l, 2])
        x_l_new = (y_l - yn)*cosA + (x_l - xn)*sinA
        y_l_new = -(y_l - yn)*sinA + (x_l - xn)*cosA
        dict_x.append(x_l_new)
        dict_y.append(y_l_new)
    print(dict_x, dict_y)
    # coords = get_polygon_coords(file.dropna)
    # transformed_coords = transform_coordinates(coords, start, stop)
    # dict_x, dict_y = zip(*transformed_coords)
    # print(dict_x, dict_y)
    # Вхождение высоты в конкретный интервал
    dict_ins = []
    for k in range(len_part - 1):
        if (H < 0 and ((H <= dict_y[k] and H >= dict_y[k+1])
                       or (H >= dict_y[k] and H <= dict_y[k+1]))):
            s = [k, k+1]
            dict_ins.append(s)
        elif (H > 0 and ((H >= dict_y[k] and H <= dict_y[k+1])
                         or (H <= dict_y[k] and H >= dict_y[k+1]))):
            s = [k, k+1]
            dict_ins.append(s)

    # Функция для определения пересечения высоты
    row_3 = file[file['STATUS'].isnull()]  # !!!!!
    for i in range(len(dict_ins)):
        x_k = dict_x[dict_ins[i][0]]
        y_k = dict_y[dict_ins[i][0]]
        x_k_1 = dict_x[dict_ins[i][1]]
        y_k_1 = dict_y[dict_ins[i][1]]
        X_H = (H - y_k)*(x_k_1 - x_k)/(y_k_1 - y_k) + x_k
        Y_H_main = X_H*cosA - H*sinA + yn
        X_H_main = X_H*sinA + H*cosA + xn
        start_row.loc[1] = [status + i, X_H_main, Y_H_main, 'new']
        start_row[['number', 'X', 'Y', 'STATUS']].to_csv(
            output_dir + '/file_' + str(i) + '.csv', index=False, sep=';')
        stop_row.to_csv(output_dir + '/file_' + str(i) + '.csv',
                        index=False, sep=';', header=None, mode='a')
        row_3.to_csv(output_dir + '/file_' + str(i) + '.csv',
                     index=False, sep=';', header=None, mode='a')
