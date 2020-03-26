from os import path
import numpy as np
import pandas as pd


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
    data = pd.read_csv(filename, sep=';', decimal=',')
    return data.astype({"x": float, "y": float})


def get_polygon_coords(dataframe, inverse=False):
    columns = ["x", "y"]
    if inverse:
        columns.reverse()
    coords = dataframe[columns].values.tolist()
    return coords


def split_polygon(coords):
    total_area = area_by_shoelace(coords)
    half_area = total_area*0.5
    triangles = [[coords[0], *coords[i:i+2]] for i in range(len(coords)-2)]
    triangle_areas = {i: area_by_shoelace(
        triangle) for i, triangle in enumerate(triangles)}
    current_area = 0
    index = -1

    while(current_area < half_area):
        index += 1
        current_area += triangle_areas.get(index)

    area = triangle_areas.get(index) - (current_area-half_area)
    splitting_point = find_splitting_point(
        triangles[index], area)
    split_line = [coords[0], splitting_point]
    new_coords = [*coords[:index+1], splitting_point, *coords[index+1:]]
    return [new_coords, split_line]


def transform_coordinates(coords, start, stop, versa=False):
    basis = dist(start, stop)
    cos, sin = np.divide(np.subtract(stop, start), basis)
    if versa == False:
        transform_matrix = np.array([[cos, sin], [-sin, cos]])
        transformed = [np.matmul(transform_matrix, np.subtract(
            point, start)).tolist() for point in coords]
        return transformed
    elif versa == True:
        transformed = [(point[0]*sin + point[1]*cos + start[0], 
                point[0]*cos - point[1]*sin + start[1]) for point in coords]
        return transformed        

      
def smooth_polygon(data, status, output_dir):
    start_row = data[data.status.eq("start")]
    stop_row = data[data.status.eq("stop")]
    start = start_row[["x", "y"]].values[0]
    stop = stop_row[["x", "y"]].values[0]
    # data for new square
    new_data = data[data.status.ne('+')]

    coords = get_polygon_coords(data)
    new_coords = get_polygon_coords(new_data)

    # difference of squares
    area_diff = area_by_shoelace(coords)-area_by_shoelace(new_coords)

    # def - для локальной СК, чтобы сразу считать
    basis = dist(start, stop)
    H = 2*area_diff/basis

    part = data.dropna()
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
    new_data = data[data.status.isnull()]
    start_row = start_row.values[0].tolist()
    stop_row = stop_row.values[0].tolist()
    for i, point in enumerate(dict_ins):
        x_k = dict_x[point[0]]
        y_k = dict_y[point[0]]
        x_k_1 = dict_x[point[1]]
        y_k_1 = dict_y[point[1]]
        x_H = (H - y_k)*(x_k_1 - x_k)/(y_k_1 - y_k) + x_k
        x_main, y_main = transform_coordinates([x_H, H], start, stop, versa = True)
        new_row = [status + i, x_main, y_main, 'new']
        df_data = [start_row, new_row, stop_row]
        df = pd.DataFrame(df_data, columns=data.columns).append(new_data)
        df.to_csv(path.join(output_dir, "file_{}.csv".format(i)),
                  index=False, sep=';')
