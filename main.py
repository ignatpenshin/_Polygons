from os import path
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

from helpers import (get_polygon_coords, read_polygon_data, smooth_polygon,
                     split_polygon, transform_coordinates, dist)

output_dir = "variants"
Path(output_dir).mkdir(parents=True, exist_ok=True)
polygon = read_polygon_data('polygon.csv')
smooth_polygon(file=polygon, status=900, output_dir=output_dir)
new_output_dir = output_dir+"_2"
smooth_polygon_file = path.join(new_output_dir, "file_31.csv")
smoothed_polygon = read_polygon_data(smooth_polygon_file)

smooth_polygon(file=smoothed_polygon, status=500, output_dir=new_output_dir)
final_polygon = path.join(new_output_dir, "file_0.csv")


dataframe = read_polygon_data(final_polygon)
coords = get_polygon_coords(dataframe)
new_polygon, split_line = split_polygon(coords)
start, stop = split_line
new_coords = transform_coordinates(new_polygon, start, stop)
start_index = new_polygon.index(start)
stop_index = new_polygon.index(stop)
# new_coords = [[0, 0], [0, 6], [3, 6], [3, 0], [3, -6], [0, -6]]
x0, y0 = new_coords[start_index]
x1, y1 = new_coords[stop_index]
x_a, y_a = new_coords[start_index-1]
x_b, y_b = new_coords[start_index+1]
x_c, y_c = new_coords[stop_index-1]
x_d, y_d = new_coords[stop_index+1]

k2 = (x_b*y0 - x0*y_b)/(x_b-x0)
a2 = (y_b - y0)/(x_b - x0)
k5 = (x_a*y0 - x0*y_a)/(x_a-x0)
a5 = (y_a - y0)/(x_a - x0)
k3 = (x_c*y_d - x_d*y_c)/(x_c-x_d)
a3 = (y_c - y_d)/(x_c - x_d)
# f1 = (1/A2 - 1/A5)
# f2 = (-(K2/A2) + K3/A3 - (24+K3)/A3 + (24+K5)/A5)
# f3 = (-(K3*y1/A3) - 12*x1 + (144 + 12*K3)/A3 +
#       (12*y1 + K3*y1)/A3 - (144 + 12*K5)/A5)

# a_1 = (-f2 + (f2**2 - 4*f1*f3)**0.5)/(2*f1)
# a_2 = (-f2 - (f2**2 - 4*f1*f3)**0.5)/(2*f1)
# print(a_1, a_2)

a_1 = -(((-a2*a3*(k5+24)+24*a2*a5 + a3*a5*k2)**2 - 48*a2*a3*(a2-a5)*(a3*(a5*x1+k5+12) -
                                                                     a5*(k3+y1+12))) - a2*a3*k5 - 24*a2*a3 + 24*a2*a5 + a3*a5*k2)**0.5/(2*a3*(a2-a5))

a_2 = -(((-a2*a3*(k5+24)+24*a2*a5 + a3*a5*k2)**2 - 48*a2*a3*(a2-a5)*(a3*(a5*x1+k5+12) -
                                                                     a5*(k3+y1+12))) + a2*a3*k5 + 24*a2*a3 - 24*a2*a5 - a3*a5*k2)**0.5/(2*a3*(a2-a5))
# print(a_1, a_2)☺

fig, ax = plt.subplots(2)
polygon = Polygon(coords)
p = PatchCollection([polygon], fc="none", ec="purple")
ax[0].add_collection(p)
x, y = zip(*new_polygon)
ax[1].set_xlim(min(x)-10,max(x)+10)
ax[1].set_ylim(min(y)-10,max(y)+10)
ax[0].plot(*zip(*split_line), "--")
ax[0].scatter(x, y, c="blue")
polygon = Polygon(new_coords)
p = PatchCollection([polygon], fc="none", ec="purple")
ax[1].add_collection(p)
x, y = zip(*new_coords)
ax[1].set_xlim(min(x)-10,max(x)+10)
ax[1].set_ylim(min(y)-10,max(y)+10)
split_line = [new_coords[start_index], new_coords[stop_index]]
ax[1].plot(*zip(*split_line), "--")
ax[1].scatter(x, y, c="blue")
plt.show()