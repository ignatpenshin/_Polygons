from pandas import read_csv
import numpy as np


def area_by_shoelace(coords):
    x, y = zip(*coords)
    return 0.5 * np.abs(np.dot(x[:-1], y[1:]) + x[-1]*y[0] -
                        np.dot(y[:-1], x[1:]) - y[-1]*x[0])


def dist(point1, point2):
    return np.linalg.norm(np.subtract(point2, point1))


def find_splitting_point(triangle, area):
    p1, p2, p3 = triangle
    k = p2[0]*(p1[1]-p3[1]) + p1[0]*(p3[1]-p2[1]) + p3[0]*(p2[1]-p1[1])
    x = (2*area*(p3[0]-p2[0])+p2[0]*k) / k
    y = (p2[1]-p3[1])*(x-p3[0])/(p2[0]-p3[0]) + p3[1]
    return [x, y]


def read_polygon_data(filename):
    file = read_csv(filename, sep=';', decimal=',')
    return file.astype({"x": float, "y": float})


def get_polygon_coords(dataframe, inverse=False):
    columns = ["x", "y"]
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


def transform_coordinates(coords, start, stop):
    basis = dist(start, stop)
    xn, yn = start
    cos, sin = np.divide(np.subtract(stop, start), basis)

    transformed = [[
        (x - xn)*cos + (y - yn)*sin, -(x - xn)*sin + (y - yn)*cos] for x, y in coords]

    return transformed


def smooth_polygon(file, status, output_dir):
    start_row = file[file.status.eq("start")]
    stop_row = file[file.status.eq("stop")]
    start = start_row[["x", "y"]].values[0]
    stop = stop_row[["x", "y"]].values[0]
    xn, yn = start
    xn_j, yn_j = stop
    # data for new square
    n_s = file[file.status.ne('+')]

    coords = get_polygon_coords(file)
    new_coords = get_polygon_coords(n_s)

    # difference of squares
    d_s = area_by_shoelace(coords)-area_by_shoelace(new_coords)

    # def - для локальной СК, чтобы сразу считать
    basis = dist(start, stop)
    cosA = (yn_j - yn)/basis
    sinA = (xn_j - xn)/basis
    H = 2*d_s/basis

    part = file.dropna()
    len_part = len(part)

    # Точки сглаживания в локальной СК.
    # По значению y_l_NEW можно понять, где входит H.
    coords = get_polygon_coords(part)
    dict_x, dict_y = zip(*transform_coordinates(coords, start, stop))

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
    row_3 = file[file.status.isnull()]  # !!!!!
    for i in range(len(dict_ins)):
        x_k = dict_x[dict_ins[i][0]]
        y_k = dict_y[dict_ins[i][0]]
        x_k_1 = dict_x[dict_ins[i][1]]
        y_k_1 = dict_y[dict_ins[i][1]]
        x_H = (H - y_k)*(x_k_1 - x_k)/(y_k_1 - y_k) + x_k
        y_H_main = x_H*cosA - H*sinA + yn
        x_H_main = x_H*sinA + H*cosA + xn
        start_row.loc[1] = [status + i, x_H_main, y_H_main, 'new']
        start_row[['number', 'x', 'y', 'status']].to_csv(
            output_dir + '/file_' + str(i) + '.csv', index=False, sep=';')
        stop_row.to_csv(output_dir + '/file_' + str(i) + '.csv',
                        index=False, sep=';', header=None, mode='a')
        row_3.to_csv(output_dir + '/file_' + str(i) + '.csv',
                     index=False, sep=';', header=None, mode='a')
