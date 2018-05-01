import os

import pdb
import pandas as pd

from spatial import run_all_clipping, plot_all_shapes, make_shapefile, plot_shape, extract_latlongs

rasters_to_clip = {
    "synoptic_temp": {"input_path": "Z:/mastergrids/Other_Global_Covariates/WorldClim_Temperature/1k/",
                      "input_pattern": "WorldClim_Max_Temp_[0-9]{2}_1k\.tif$",
                      "unit": "month"},
    # "baseline_incidence": {"input_path": "Z:/cubes/Pf_results/MODEL_43/output_rasters/incidence/ALL/",
    #                             "input_pattern": "MODEL43\.2015\.[0-9]{1}\.inc\.rate\.PR\.ALL\.tif$"}
}

main_path = "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/projects/Mozambique/"
out_path = os.path.join(main_path, "incidence_calibration/")
shapefile_name = "bbox"
data_path = "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/data/Mozambique/Magude"
grid_path = os.path.join(main_path, "gridded_simulation_input/grid_lookup_friction.csv")
concave_alpha = 10

print("loading and cleaning point data")
grid_data = pd.read_csv(grid_path)
grid_data.rename(columns={"mid_x": "lng", "mid_y": "lat"}, inplace=True)

for input_name, input_vals in rasters_to_clip.items():

    print("extracting for " + input_name)

    # plot_all_shapes(hh_data, alphas=[1, 5, 10])
    run_all_clipping(input_vals["input_path"], out_path, shapefile_name, grid_data,
                     raster_pattern=input_vals["input_pattern"], overwrite=True,
                     raster_folder=input_name, alpha=concave_alpha, out_name="{name}_all".format(name=input_name),
                     write_shp=True, unit=input_vals["unit"])

    clip_catchments = False
    if clip_catchments:

        for  hf in pd.unique(grid_data['catchment']):
            run_all_clipping(os.path.join(out_path, "rasters", input_name), out_path, shapefile_name,
                             grid_data.query('catchment==@hf'), raster_pattern="{name}_all.*tif$".format(name=input_name),
                             overwrite=True, raster_folder=input_name, alpha=concave_alpha,
                             out_name="{name}_{hf_name}_{shp_type}".format(name=input_name, hf_name=hf, shp_type=shapefile_name),
                             write_shp=False, crop=False)

