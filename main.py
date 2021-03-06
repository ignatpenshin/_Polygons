from helpers import (get_polygon_coords, read_polygon_data, smooth_polygon, area_by_shoelace,
                     split_polygon, transform_coordinates, dist, road_project)
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from pathlib import Path
from os import path
import pandas as pd 

ROAD_WIDTH = 12
output_dir = "variants"
Path(output_dir).mkdir(parents=True, exist_ok=True)
polygon = read_polygon_data('polygon.csv')
smooth_polygon(data=polygon, status=900, output_dir=output_dir)
new_output_dir = output_dir+"_2"
smooth_polygon_file = path.join(new_output_dir, "file_31.csv")
smoothed_polygon = read_polygon_data(smooth_polygon_file)



smooth_polygon(data=smoothed_polygon, status=500, output_dir=new_output_dir)
final_polygon = path.join(new_output_dir, "file_0.csv")

road = road_project(final_polygon, 12,0)


print(pd.DataFrame({'x':r[0], 'y':r[1]} for r in road))

dataframe = read_polygon_data(final_polygon)
coords = get_polygon_coords(dataframe, inverse=True)
new_polygon, split_line = split_polygon(coords, 0)
start, stop = split_line
transformed_coords = transform_coordinates(new_polygon, start, stop)
start_index = new_polygon.index(start)
stop_index = new_polygon.index(stop)


x0, y0 = transformed_coords[start_index]
x_1, y_1 = transformed_coords[stop_index]
x_a, y_a = transformed_coords[start_index-1]
x_b, y_b = transformed_coords[start_index+1]
x_c, y_c = transformed_coords[stop_index-1]
x_d, y_d = transformed_coords[stop_index+1]

h = (((-2*x_1*y_a*y_b*y_c+2*x_1*y_a*y_b*y_d-24*x_a*y_b*y_c+24*x_a*y_b*y_d+2*x_c*y_a*y_b*y_d+24*x_c*y_a*y_b-2*x_d*y_a*y_b*y_c-24*x_d*y_a*y_b)**2
      - 4*(ROAD_WIDTH*x_1*y_a*y_b*y_c-ROAD_WIDTH*x_1*y_a*y_b*y_d+144*x_a*y_b*y_c-144*x_a*y_b*y_d-ROAD_WIDTH*x_c*y_a*y_b*y_d-144*x_c*y_a*y_b+ROAD_WIDTH*x_d*y_a*y_b*y_c
           + 144*x_d*y_a*y_b)*(x_a * y_b*y_c-x_a*y_b*y_d+x_b*y_a*y_c-x_b*y_a*y_d-2*x_c*y_a*y_b+2*x_d*y_a*y_b))**0.5
     + 2*x_1*y_a*y_b*y_c-2*x_1*y_a*y_b*y_d+24*x_a*y_b*y_c-24*x_a*y_b*y_d-2*x_c*y_a*y_b*y_d-24*x_c*y_a*y_b+2*x_d*y_a*y_b*y_c+24*x_d*y_a*y_b)/(2*(x_a*y_b*y_c
                                                                                                                                                - x_a*y_b*y_d+x_b*y_a*y_c-x_b*y_a*y_d-2*x_c*y_a*y_b+2*x_d*y_a*y_b))
h2 = (-((-2*x_1*y_a*y_b*y_c+2*x_1*y_a*y_b*y_d-24*x_a*y_b*y_c+24*x_a*y_b*y_d+2*x_c*y_a*y_b*y_d+24*x_c*y_a*y_b-2*x_d*y_a*y_b*y_c-24*x_d*y_a*y_b)**2
        - 4*(ROAD_WIDTH*x_1*y_a*y_b*y_c-ROAD_WIDTH*x_1*y_a*y_b*y_d+144*x_a*y_b*y_c-144*x_a*y_b*y_d-ROAD_WIDTH*x_c*y_a*y_b*y_d-144*x_c*y_a*y_b+ROAD_WIDTH*x_d*y_a*y_b*y_c
             + 144*x_d*y_a*y_b)*(x_a * y_b*y_c-x_a*y_b*y_d+x_b*y_a*y_c-x_b*y_a*y_d-2*x_c*y_a*y_b+2*x_d*y_a*y_b))**0.5
      + 2*x_1*y_a*y_b*y_c-2*x_1*y_a*y_b*y_d+24*x_a*y_b*y_c-24*x_a*y_b*y_d-2*x_c*y_a*y_b*y_d-24*x_c*y_a*y_b+2*x_d*y_a*y_b*y_c+24*x_d*y_a*y_b)/(2*(x_a*y_b*y_c
                                                                                                                                                 - x_a*y_b*y_d+x_b*y_a*y_c-x_b*y_a*y_d-2*x_c*y_a*y_b+2*x_d*y_a*y_b))

x_2 = h*x_b/y_b
x_3 = (h-y_d)*(x_c-x_d)/(y_c-y_d)+x_d
x_4 = (h-ROAD_WIDTH-y_d)*(x_c-x_d)/(y_c-y_d)+x_d
x_5 = (h-ROAD_WIDTH)*x_a/y_a

mid_top = (x_3+x_2)/2
top_polygon = [[mid_top, h], [x_2, h], *
               transformed_coords[start_index+1:stop_index], [x_3, h]]
new_polygon_top, split_line_top = split_polygon(top_polygon, 0)
new_polygon_top_1 = new_polygon_top[:new_polygon_top.index(
    split_line_top[1])+1]
new_polygon_top_2 = [*new_polygon_top[new_polygon_top.index(
    split_line_top[1]):], new_polygon_top[0]]

mid_bottom = (x_5+x_4)/2
bottom_polygon = [[mid_bottom, h-ROAD_WIDTH], [x_4, h-ROAD_WIDTH], *
                  transformed_coords[stop_index+1:], [x_5, h-ROAD_WIDTH]]
new_polygon_bottom, split_line_bottom = split_polygon(bottom_polygon,0)
new_polygon_bottom_1 = new_polygon_bottom[:new_polygon_bottom.index(
    split_line_bottom[1])+1]
new_polygon_bottom_2 = [*new_polygon_bottom[new_polygon_bottom.index(
    split_line_bottom[1]):], new_polygon_bottom[0]]

road_top = [[x_3, h], transformed_coords[stop_index],
            transformed_coords[start_index], [x_2, h]]
road_bottom = [transformed_coords[stop_index], [x_4, h-ROAD_WIDTH],
               [x_5, h-ROAD_WIDTH], transformed_coords[start_index]]

result_polygons = [new_polygon_top_1, new_polygon_top_2,
                   new_polygon_bottom_1, new_polygon_bottom_1, road_top, road_bottom]

result_area = sum(
    map(lambda polygon: area_by_shoelace(polygon), result_polygons))
# print("initial total area: {}\nresult total area: {}".format(area_by_shoelace(coords), result_area)
#       )

x_6 = split_line_top[1]
x_7 = split_line_bottom[1]

cut_points = [top_polygon[0], bottom_polygon[0],
              split_line_bottom[1], split_line_top[1], [x_2, h], [x_3, h], [x_4, h-ROAD_WIDTH], [x_5, h-ROAD_WIDTH]]

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
polygon = Polygon(transformed_coords)
p = PatchCollection([polygon], fc="none", ec="grey")
ax[1].add_collection(p)
x, y = zip(*transformed_coords)
ax[1].axis("equal")
split_line = [transformed_coords[start_index], transformed_coords[stop_index]]
ax[1].plot(*zip(*split_line), "--")
ax[1].scatter(x, y, c="blue")
polygon = Polygon(new_polygon_top_1)
p = PatchCollection([polygon], fc="none", ec="green")
ax[1].add_collection(p)
polygon = Polygon(new_polygon_top_2)
p = PatchCollection([polygon], fc="none", ec="red")
ax[1].add_collection(p)
polygon = Polygon(new_polygon_bottom_1)
p = PatchCollection([polygon], fc="none", ec="cyan")
ax[1].add_collection(p)
polygon = Polygon(new_polygon_bottom_2)
p = PatchCollection([polygon], fc="none", ec="magenta")
ax[1].add_collection(p)
x, y = zip(*cut_points)
ax[1].scatter(x, y, c="orange")
plt.show()
