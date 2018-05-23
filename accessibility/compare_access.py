## compare Dan Weiss' 2015 "friction surface" with Caitlin Bever's health seeking function

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import rasterio
import seaborn as sns
import shapely.geometry
import pdb

from spatial import  make_shapefile, extract_latlongs

main_dir = "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/projects/Mozambique"
weiss_dir =os.path.join(main_dir,"incidence_calibration/rasters/accessibility_to_cities/accessibility_to_cities_all.year.2015.tif" )
weiss_dir =os.path.join(main_dir,"incidence_calibration/rasters/magude_access/study.area.accessibility.tif" )
bever_dir = os.path.join(main_dir, "gridded_simulation_input")
bever_var = 'cov_newclin_youth'

print("loading bever surface")
bever_surface = pd.read_csv(os.path.join(bever_dir, "grid_all_healthseek_events.csv"))
bever_lookup = pd.read_csv(os.path.join(bever_dir, "grid_lookup.csv"))
healthseek_compare = pd.merge(bever_surface[["grid_cell", bever_var]], bever_lookup, how="outer")

sim_shp = make_shapefile(healthseek_compare.copy(), type='point',
                            to_crs={'init' :'epsg:4326'},
                            lat_name='mid_y', lon_name='mid_x')
shapes = [shapely.geometry.mapping(g) for g in sim_shp['geometry']]

healthseek_compare['weiss_access'] = extract_latlongs(weiss_dir, shapes)
healthseek_compare = healthseek_compare.query('weiss_access>0')

print("plotting")
fig= plt.figure(figsize=(11,10))

bever_ax = plt.subplot(221)
weiss_ax = plt.subplot(222)
cor_ax = plt.subplot(223)
full_weiss_ax = plt.subplot(224)

print("plotting bever")
# bever_ax = axes[0]
ab = bever_ax.scatter(healthseek_compare['mid_x'], healthseek_compare['mid_y'],
                      c=healthseek_compare['cov_newclin_youth'], marker=',')
bever_ax.set_title("Bever: P(Health-Seeking)")
bever_ax.xaxis.set_visible(False)
bever_ax.yaxis.set_visible(False)
plt.colorbar(ab, ax=bever_ax)

print("plotting weiss")
# weiss_ax = axes[1]
aw = weiss_ax.scatter(healthseek_compare['mid_x'], healthseek_compare['mid_y'],
                      c=healthseek_compare['weiss_access'], vmin=0, marker=',', cmap='viridis_r')
plt.colorbar(aw, ax=weiss_ax)
weiss_ax.set_title("Weiss: Accessibility (min)")
weiss_ax.xaxis.set_visible(False)
weiss_ax.yaxis.set_visible(False)

print("plotting correlation")
for_corr = healthseek_compare.query('weiss_access>0')
# cor_ax = axes[2]
corr = cor_ax.scatter(for_corr['cov_newclin_youth'], for_corr['weiss_access'],
                      c=for_corr['cov_newclin_youth'], vmin=0)
plt.colorbar(ab, ax=cor_ax)
cor_ax.set_title("Correlation between variables")
cor_ax.set_xlabel("P(Health-Seeking)")
cor_ax.set_ylabel("Weiss Accessibility (min)")

print("plotting full Weiss heatmap")
weiss_surface = rasterio.open(weiss_dir).read(1)
masked_array = np.ma.array (weiss_surface, mask=weiss_surface<0)
cmap = mpl.cm.viridis_r
cmap.set_bad('white',1.)
full_weiss = full_weiss_ax.imshow(masked_array, vmin=0, cmap=cmap)
plt.colorbar(full_weiss, ax=full_weiss_ax)
full_weiss_ax.set_title("Weiss Accessibility, Magude")
full_weiss_ax.xaxis.set_visible(False)
full_weiss_ax.yaxis.set_visible(False)

fig.tight_layout()

plt.savefig("hf_accessibility_compare.png")


plt.show()
