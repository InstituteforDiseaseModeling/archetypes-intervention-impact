import pandas as pd

from covariate_prep.spatial import run_all_clipping, plot_all_shapes

# setup
tsi_path = "Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/1k/Synoptic-Months"
tsi_pattern = ".*Mean.*tif$"

# map_covariate_path = "Z:/mastergrids/MODIS_Global/MCD12Q1_Annual_Landcover/500m_Raw"
# covariate_pattern = re.compile(".*tif$")
out_path = "Z:/amelia/larval_habitat/karen/"
shapefile_name = "point"
data_path = "Q:/Malaria/Myanmar/Karen/metf_data/prevalence_data.csv"

# make point dataframe with three columns: id, latitude, longitude
print("loading and cleaning point data")
data = pd.read_csv(data_path, encoding='mac_roman')

plot_all_shapes(data)

run_all_clipping(tsi_path, out_path, shapefile_name, data,
                 raster_pattern=tsi_pattern, overwrite=False,
                 raster_folder="tsi")
