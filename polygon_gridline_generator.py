import geopandas as gpd
from shapely.geometry import Polygon, Point
from shapely.affinity import rotate
import numpy as np
import webbrowser
import os
import math
import pandas as pd
import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


def calculate_tilt(p1, p2):
    return math.degrees(math.atan2(p2[1]-p1[1], p2[0]-p1[0]))

def rotate_polygon(poly, angle, origin):
    return rotate(poly, angle, origin=origin, use_radians=False)

def generate_grid_for_polygon(polygon, poly_name, start_labels, first_row_inc, color_dict, crs):
    height = 33 * 16 / 3.28  
    width  = 33 * 10 / 3.28  

    pts = list(polygon.exterior.coords)
    angle = calculate_tilt(pts[0], pts[1])
    minx, miny = pts[0][0], pts[0][1]

    grid_cells, placemark_points, labels = [], [], []
    x_coords = np.arange(minx, minx + 5000, width)
    y_coords = np.arange(miny, miny + 5000, height)
    if not start_labels: 
        start_labels = [1]
    
    for j, y in enumerate(y_coords):
        if j >= len(start_labels): 
            break
        label = start_labels[j]
        for x in x_coords:
            cell = Polygon([(x, y), (x+width, y), (x+width, y+height), (x, y+height)])
            if angle != 0:
                cell = rotate_polygon(cell, angle, Point(minx, miny))
            if polygon.intersects(cell):
                clip = polygon.intersection(cell)
                if not clip.is_empty:
                    grid_cells.append(clip)
                    placemark_points.append(Point(clip.centroid.x, clip.centroid.y))
                    labels.append(label)
            label += 1 if (j % 2 == 0) == first_row_inc else -1

    grid_gdf = gpd.GeoDataFrame({'name': labels, 'geometry': grid_cells}, crs=crs)
    placemarks_gdf = gpd.GeoDataFrame({'name': labels, 'geometry': placemark_points}, crs=crs)
    grid_gdf['polygon_name'] = poly_name
    placemarks_gdf['polygon_name'] = poly_name
    grid_gdf['color_dict'] = [color_dict]*len(grid_gdf)
    return grid_gdf, placemarks_gdf

def save_and_open_kml(all_grid, all_places, output_kml, settings):
    if not all_grid.crs.is_geographic:
        all_grid['area_m2'] = all_grid.geometry.area
    else:
        proj = all_grid.to_crs(all_grid.estimate_utm_crs())
        all_grid['area_m2'] = proj.geometry.area

    color_map = {
    "green":   "8000ff00",
    "red":     "800000ff",
    "blue":    "80ff0000",
    "yellow":  "8000ffff",
    "purple":  "80800080",
    "orange":  "8000a5ff",
    "cyan":    "80ffff00",
    "magenta": "80ff00ff",
    "brown":   "802a2aa5",
    "gray":    "80808080" 
}

    if not all_grid.crs.is_geographic:
        all_grid = all_grid.to_crs(epsg=4326)
        all_places = all_places.to_crs(epsg=4326)
    
    with open(output_kml, 'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n')
        for idx, row in all_grid.iterrows():
            label = row["name"]
            area_acres = row['area_m2'] / 4046.86
            acre = int(area_acres)
            gunta = round((area_acres - acre) * 40)
            if gunta == 40:
                acre += 1
                gunta = 0
            area_txt = f" ({acre} - {gunta})"
            poly_name = row.get('polygon_name','')
            sett = settings.get(poly_name, {})
            cdict = sett.get('color_dict', {})
            has_color, cell_color, width_line = False, "800000ff", 2
            for col, labs in cdict.items():
                if label in labs:
                    cell_color = color_map.get(col, "800000ff")
                    has_color = True
                    break
            poly_fill = "1" if has_color else "0"
            f.write(f'<Placemark><name>{label}{area_txt}</name>\n')
            f.write(f'<Style><PolyStyle><fill>{poly_fill}</fill><color>{cell_color}</color></PolyStyle>')
            f.write(f'<LineStyle><color>800000ff</color><width>{width_line}</width></LineStyle></Style>\n')
            f.write('<Polygon><outerBoundaryIs><LinearRing><coordinates>\n')
            for coord in row.geometry.exterior.coords:
                f.write(f'{coord[0]},{coord[1]},0\n')
            f.write('</coordinates></LinearRing></outerBoundaryIs></Polygon>\n</Placemark>\n')
        for idx, row in all_places.iterrows():
            area_acres = all_grid.iloc[idx]['area_m2'] / 4046.86
            acre = int(area_acres)
            gunta = round((area_acres - acre) * 40)
            if gunta == 40:
                acre += 1
                gunta = 0
            area_txt = f" ({acre} - {gunta})"
            if (acre == 4 and gunta == 0):
                area_txt = ""
            f.write(f'<Placemark><name>{row["name"]}{area_txt}</name>\n')
            f.write('<Style><IconStyle><scale>0</scale></IconStyle>')
            f.write('<LabelStyle><scale>1</scale></LabelStyle></Style>\n')
            f.write('<Point><coordinates>\n')
            f.write(f'{row.geometry.x},{row.geometry.y},0\n')
            f.write('</coordinates></Point>\n</Placemark>\n')
        f.write('</Document>\n</kml>\n')
    webbrowser.open(f"file://{os.path.realpath(output_kml)}")

def generate_grids_for_multiple_polygons(kml_path, settings, output_kml="combined_grid_output.kml", extensions=None):
    gdf = gpd.read_file(kml_path)
    if gdf.empty: 
        raise ValueError("No polygons found.")
    if gdf.crs.is_geographic: 
        gdf = gdf.to_crs(gdf.estimate_utm_crs())
    
    grid_list, place_list = [], []
    for idx, row in gdf.iterrows():
        poly = row.geometry
        poly_name = row.get('Name', f"Polygon_{idx}")
        sett = settings.get(poly_name, {})
        start_labels = sett.get('start_labels', [1])
        first_row_inc = sett.get('first_row_increment', True)
        color_dict = sett.get('color_dict', None)
        grid_gdf, placemarks_gdf = generate_grid_for_polygon(poly, poly_name, start_labels, first_row_inc, color_dict, gdf.crs)
        grid_list.append(grid_gdf)
        place_list.append(placemarks_gdf)
        print(f"Grid generated for polygon: {poly_name}")
    combined_grid = gpd.GeoDataFrame(pd.concat(grid_list, ignore_index=True), crs=gdf.crs)
    combined_places = gpd.GeoDataFrame(pd.concat(place_list, ignore_index=True), crs=gdf.crs)
    
    if extensions:
        combined_grid = apply_multiple_extensions(combined_grid, extensions)
    save_and_open_kml(combined_grid, combined_places, output_kml, settings)

def find_extreme_edge(poly, side):
    coords = list(poly.exterior.coords)
    if coords[0] == coords[-1]:
        coords = coords[:-1]
    n = len(coords)
    extreme_val, edge_index = None, None
    if side in ("bottom", "top"):
        for i in range(n):
            j = (i+1) % n
            avg = (coords[i][1] + coords[j][1]) / 2.0
            if extreme_val is None or (side=="bottom" and avg < extreme_val) or (side=="top" and avg > extreme_val):
                extreme_val = avg
                edge_index = i
    elif side in ("left", "right"):
        for i in range(n):
            j = (i+1) % n
            avg = (coords[i][0] + coords[j][0]) / 2.0
            if extreme_val is None or (side=="left" and avg < extreme_val) or (side=="right" and avg > extreme_val):
                extreme_val = avg
                edge_index = i
    return edge_index, coords

def extend_polygon_with_edge(source_poly, target_poly, direction):
    mapping = {
        "bottom_to_top": ("bottom", "top"),
        "top_to_bottom": ("top", "bottom"),
        "left_to_right": ("left", "right"),
        "right_to_left": ("right", "left")
    }
    if direction not in mapping:
        raise ValueError("Invalid direction")
    s_side, t_side = mapping[direction]
    s_index, s_coords = find_extreme_edge(source_poly, s_side)
    t_index, t_coords = find_extreme_edge(target_poly, t_side)
    # Get target edge vertices
    t1 = t_coords[t_index]
    t2 = t_coords[(t_index+1) % len(t_coords)]
    # Ensure insertion order (adjust based on axis)
    if s_side in ("bottom", "top"):
        s_pt1, s_pt2 = s_coords[s_index], s_coords[(s_index+1)%len(s_coords)]
        if s_pt1[0] < s_pt2[0] and t1[0] > t2[0]:
            t1, t2 = t2, t1
        elif s_pt1[0] > s_pt2[0] and t1[0] < t2[0]:
            t1, t2 = t2, t1
    else:
        s_pt1, s_pt2 = s_coords[s_index], s_coords[(s_index+1)%len(s_coords)]
        if s_pt1[1] < s_pt2[1] and t1[1] > t2[1]:
            t1, t2 = t2, t1
        elif s_pt1[1] > s_pt2[1] and t1[1] < t2[1]:
            t1, t2 = t2, t1

    new_coords = s_coords[:s_index+1] + [t1, t2] + s_coords[s_index+1:]
    if new_coords[0] != new_coords[-1]:
        new_coords.append(new_coords[0])
    return Polygon(new_coords)

def extend_polygon_in_gdf(gdf, source_label, target_label, direction):
    src_rows = gdf[gdf['name'] == source_label]
    tgt_rows = gdf[gdf['name'] == target_label]
    if src_rows.empty or tgt_rows.empty:
        print("Source or target polygon not found.")
        return gdf
    src_idx = src_rows.index[0]
    new_poly = extend_polygon_with_edge(gdf.at[src_idx, 'geometry'], tgt_rows.iloc[0].geometry, direction)
    gdf.at[src_idx, 'geometry'] = new_poly
    return gdf

def apply_multiple_extensions(gdf, extension_instructions):
    for ext in extension_instructions:
        gdf = extend_polygon_in_gdf(gdf, ext["source_label"], ext["target_label"], ext["direction"])
    return gdf


KML_FILE = "Polygon.kml"
OUTPUT_KML = "combined_grid_output.kml"

# Define your settings and extension instructions, e.g.:
settings_by_polygon = {
    "A": {"start_labels": [479,478,428,427,386,386,347,343,306,309,272],
          "first_row_increment": True,
          "color_dict": {"green": [479,478,428]}
          },
    "B": {"start_labels": [427,386,386,347,343,306,309,272],
          "first_row_increment": False,
          "color_dict": {"red": [427,386,386]}
          },
}

extensions = [
    # {"source_label": 479, "target_label": 697, "direction": "left_to_right"},
    # {"source_label": 478, "target_label": 642, "direction": "left_to_right"},
    # {"source_label": 428, "target_label": 613, "direction": "left_to_right"},
]

class KMLFileEventHandler(FileSystemEventHandler):
    def __init__(self, kml_path, settings, output_kml, extensions=None):
        if extensions is None:
            extensions = []  
        self.kml_path = os.path.abspath(kml_path)
        self.settings = settings
        self.output_kml = output_kml
        self.extensions = extensions

    def on_modified(self, event):
        # Check if the modified file is our KML file
        if os.path.abspath(event.src_path) == self.kml_path:
            print(f"{self.kml_path} modified. Regenerating grids...")
            try:
                generate_grids_for_multiple_polygons(self.kml_path, self.settings, self.output_kml, self.extensions)
                print("Grids regenerated successfully.")
            except Exception as e:
                print(f"Error generating grids: {e}")
        else:
            print(f"Ignored change: {event.src_path}")  # Debugging: when it's not the correct file

if __name__ == "__main__":
    KML_FILE = "Polygon.kml"  # Define the correct path here

    if not os.path.exists(KML_FILE):
        print(f"Error: {KML_FILE} does not exist in the folder.")
    else:
        event_handler = KMLFileEventHandler(KML_FILE, settings_by_polygon, OUTPUT_KML)
        observer = Observer()
        
        # Watch the directory containing the KML file
        observer.schedule(event_handler, path=os.path.dirname(os.path.abspath(KML_FILE)), recursive=False)
        observer.start()
        
        print(f"Watching for changes in {KML_FILE} ...")
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            observer.stop()
        observer.join()
