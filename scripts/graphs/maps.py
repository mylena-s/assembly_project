# import folium
import json
import geopandas as gpd
import contextily as ctx
import matplotlib.pyplot as plt
from shapely.geometry import box
from matplotlib.gridspec import GridSpec
# Load GeoJSON files
with open('/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/maps/iguazu.geojson') as f:
    geojson_data = json.load(f)
with open('/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/maps/uruguay.geojson') as f:
    geojson_data2 = json.load(f)

# Filter features by attributes
iguazu = [
    feature for feature in geojson_data['features']
    if feature['properties'].get('HYBAS_ID') == 6050793620
]

iguazu = {
    "type": "FeatureCollection",
    "features": iguazu
}

uruguay = [
    feature for feature in geojson_data2['features']
    if feature['properties'].get('SUB_AREA') == 266133.3
]

uruguay = {
    "type": "FeatureCollection",
    "features": uruguay
}

# Convert GeoJSON dicts to GeoDataFrames (assumed CRS EPSG:4326)
iguazu_gdf = gpd.GeoDataFrame.from_features(iguazu['features'], crs="EPSG:4326")
uruguay_gdf = gpd.GeoDataFrame.from_features(uruguay['features'], crs="EPSG:4326")

# Define bounding box for South America in EPSG:4326
# Define square bounding box in degrees
south_america_bounds = box(-82, -40, -32, 10)  # 50x50 degrees
bbox_gdf = gpd.GeoDataFrame({'geometry': [south_america_bounds]}, crs="EPSG:4326")

# Project to EPSG:3857
bbox_gdf = bbox_gdf.to_crs(epsg=3857)
bbox_3857 = bbox_gdf.geometry.iloc[0]

# Check actual width and height
minx, miny, maxx, maxy = bbox_3857.bounds
width = maxx - minx
height = maxy - miny
print(f"Projected size: {width/1000:.1f} km x {height/1000:.1f} km")

# Plot
fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))  # square figure

iguazu_gdf = iguazu_gdf.to_crs(epsg=3857)
uruguay_gdf = uruguay_gdf.to_crs(epsg=3857)

iguazu_gdf.plot(ax=ax1, facecolor='#bb3e03', edgecolor='#bb3e03', alpha=0.8)
uruguay_gdf.plot(ax=ax1, facecolor='#ee9b00', edgecolor='#ee9b00', alpha=0.8)
# Set limits and square aspect
ax1.set_xlim(minx, maxx)
ax1.set_ylim(miny, maxy)
ax1.set_aspect('equal')  # enforce square axis
ax1.axis('off')
plt.tight_layout()# Add basemap (now will load correctly)

attribution = ""
ctx.add_basemap(ax1, source=ctx.providers.Esri.WorldPhysical, zoom=6, attribution=attribution)
plt.savefig("southam_map.svg")


################# iguazú falls

from shapely.geometry import box, Point
 
# Iguazú Falls center point (lat, lon)
garganta = gpd.GeoDataFrame(
    {'geometry': [Point(-54.436865265513546, -25.695578121699814)]},
    crs='EPSG:4326'
).to_crs(epsg=3857)

# Create GeoDataFrame for center point
center_gdf = gpd.GeoDataFrame(
    {'geometry': [Point(-54.4331745460007, -25.69603251056214)]},
    crs='EPSG:4326'
).to_crs(epsg=3857)

# Convert to Web Mercator for plotting
center_gdf = center_gdf.to_crs(epsg=3857)
center_point = center_gdf.geometry.iloc[0]

# Define bounding box with 1250m buffer
buffer_meters = 1300
bbox_iguazu = center_point.buffer(buffer_meters).envelope

# Define points in (lon, lat) order
CBFH00195 = gpd.GeoDataFrame(
    {'geometry': [Point(-54.432864, -25.692998)]},
    crs='EPSG:4326'
).to_crs(epsg=3857)

Paratype = gpd.GeoDataFrame(
    {'geometry': [Point(-54.42536, -25.70302)]},
    crs='EPSG:4326'
).to_crs(epsg=3857)



# Plot setup

# Plot points
fig, ax2 = plt.subplots(1, 1, figsize=(8, 8))

CBFH00195.plot(ax=ax2, color='#bb3e03', marker='^', label='CBFH00195',markersize=500)
Paratype.plot(ax=ax2, color='blue', marker='^', label='Paratype',markersize=500)

# Set axis limits to bounding box
ax2.set_xlim(bbox_iguazu.bounds[0], bbox_iguazu.bounds[2])
ax2.set_ylim(bbox_iguazu.bounds[1], bbox_iguazu.bounds[3])

# Add basemap with zoom appropriate for ~50 km radius
basemap = "https://tile.openstreetmap.bzh/ca/{z}/{x}/{y}.png"
attribution = ""

ctx.add_basemap(
    ax2,
    source=basemap,
    attribution=attribution,
    zoom=19)
plt.axis('off')
plt.tight_layout()
plt.savefig("Garganta_map.svg")


##### for interactive map
# # Create map centered on South America
# m = folium.Map(location=[-15, -60], zoom_start=3, tiles=None, max_zoom=8)

# # Add Esri World Physical Map tiles
# folium.TileLayer(
#     tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Physical_Map/MapServer/tile/{z}/{y}/{x}',
#     attr='Tiles © Esri — Source: US National Park Service',
#     name='Esri World Physical Map',
#     max_zoom=8,
#     overlay=False,
#     control=True
# ).add_to(m)

# # Add filtered GeoJSON layer with styling
# folium.GeoJson(
#     iguazu,
#     name='Filtered Iguazu Basin',
#     style_function=lambda feature: {
#         'fillColor': '#bb3e03',
#         'color': '#bb3e03',
#         'weight': 2,
#         'fillOpacity': 0.8,
#     }
# ).add_to(m)

# # Add filtered GeoJSON layer with styling
# folium.GeoJson(
#     uruguay,
#     name='Filtered Iguazu Basin',
#     style_function=lambda feature: {
#         'fillColor': '#ee9b00',
#         'color': '#ee9b00',
#         'weight': 2,
#         'fillOpacity': 0.8,
#     }
# ).add_to(m)

# # Add layer control
# folium.LayerControl().add_to(m)

# # Save map
# m.save('southamerica_filtered_iguazu_basin.html')
