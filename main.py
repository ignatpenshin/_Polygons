from helpers import (get_polygon_coords, read_polygon_data, smooth_polygon, area_by_shoelace,
                     split_polygon, transform_coordinates, dist)
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from pathlib import Path
from os import path


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

x0, y0 = new_coords[start_index]
x_1, y_1 = new_coords[stop_index]
x_a, y_a = new_coords[start_index-1]
x_b, y_b = new_coords[start_index+1]
x_c, y_c = new_coords[stop_index-1]
x_d, y_d = new_coords[stop_index+1]

h = (((-2*x_1*y_a*y_b*y_c+2*x_1*y_a*y_b*y_d-24*x_a*y_b*y_c+24*x_a*y_b*y_d+2*x_c*y_a*y_b*y_d+24*x_c*y_a*y_b-2*x_d*y_a*y_b*y_c-24*x_d*y_a*y_b)**2
      - 4*(12*x_1*y_a*y_b*y_c-12*x_1*y_a*y_b*y_d+144*x_a*y_b*y_c-144*x_a*y_b*y_d-12*x_c*y_a*y_b*y_d-144*x_c*y_a*y_b+12*x_d*y_a*y_b*y_c
           + 144*x_d*y_a*y_b)*(x_a * y_b*y_c-x_a*y_b*y_d+x_b*y_a*y_c-x_b*y_a*y_d-2*x_c*y_a*y_b+2*x_d*y_a*y_b))**0.5
     + 2*x_1*y_a*y_b*y_c-2*x_1*y_a*y_b*y_d+24*x_a*y_b*y_c-24*x_a*y_b*y_d-2*x_c*y_a*y_b*y_d-24*x_c*y_a*y_b+2*x_d*y_a*y_b*y_c+24*x_d*y_a*y_b)/(2*(x_a*y_b*y_c
                                                                                                                                                - x_a*y_b*y_d+x_b*y_a*y_c-x_b*y_a*y_d-2*x_c*y_a*y_b+2*x_d*y_a*y_b))


x_2 = h*x_b/y_b
x_3 = (h-y_d)*(x_c-x_d)/(y_c-y_d)+x_d
x_5 = (h-12)*x_a/y_a
x_4 = (h-12-y_d)*(x_c-x_d)/(y_c-y_d)+x_d
mid_top = (x_3+x_2)/2
mid_bottom = (x_5+x_4)/2
top_polygon = [[mid_top, h], [x_2, h], *
               new_coords[start_index+1:stop_index], [x_3, h]]
new_polygon_top, split_line_top = split_polygon(top_polygon)
new_polygon_top_1 = new_polygon_top[:new_polygon_top.index(
    split_line_top[1])+1]
new_polygon_top_2 = [*new_polygon_top[new_polygon_top.index(
    split_line_top[1]):], new_polygon_top[0]]
bottom_polygon = [[mid_bottom, h-12], [x_4, h-12], *
                  new_coords[stop_index+1:], [x_5, h-12]]
new_polygon_bottom, split_line_bottom = split_polygon(bottom_polygon)
new_polygon_bottom_1 = new_polygon_bottom[:new_polygon_bottom.index(
    split_line_bottom[1])+1]
new_polygon_bottom_2 = [*new_polygon_bottom[new_polygon_bottom.index(
    split_line_bottom[1]):], new_polygon_bottom[0]]
print(area_by_shoelace(new_polygon_top_1), area_by_shoelace(new_polygon_top_2),
      area_by_shoelace(new_polygon_bottom_1), area_by_shoelace(new_polygon_bottom_2))

fig, ax = plt.subplots(2)
polygon = Polygon(new_polygon)
p = PatchCollection([polygon], fc="none", ec="grey")
ax[0].add_collection(p)
x, y = zip(*new_polygon)
ax[0].axis("equal")
ax[0].set_xlim(min(x)-10, max(x)+10)
ax[0].set_ylim(min(y)-10, max(y)+10)
ax[0].plot(*zip(*split_line), "--")
ax[1].plot(*zip(*split_line_top), "--")
ax[1].plot(*zip(*split_line_bottom), "--")
ax[0].scatter(x, y, c="blue")
polygon = Polygon(new_coords)
p = PatchCollection([polygon], fc="none", ec="grey")
ax[1].add_collection(p)
x, y = zip(*new_coords)
ax[1].axis("equal")
split_line = [new_coords[start_index], new_coords[stop_index]]
ax[1].plot(*zip(*split_line), "--")
ax[1].scatter(x, y, c="blue")
test = [[x_3, h], new_coords[stop_index], new_coords[start_index], [x_2, h]]
# print(area_by_shoelace(test))
polygon = Polygon(new_polygon_top_1)
p = PatchCollection([polygon], fc="none", ec="green")
ax[1].add_collection(p)
test = [new_coords[stop_index], [x_4, h-12],
        [x_5, h-12], new_coords[start_index]]
# print(area_by_shoelace(test))
polygon = Polygon(new_polygon_top_2)
p = PatchCollection([polygon], fc="none", ec="red")
ax[1].add_collection(p)
polygon = Polygon(new_polygon_bottom_1)
p = PatchCollection([polygon], fc="none", ec="cyan")
ax[1].add_collection(p)
polygon = Polygon(new_polygon_bottom_2)
p = PatchCollection([polygon], fc="none", ec="magenta")
ax[1].add_collection(p)
plt.show()
