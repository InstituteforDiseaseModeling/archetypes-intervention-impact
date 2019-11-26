###############################################################################################################
## generate_era5_climate.py
## Amelia Bertozzi-Villa
## October 2019
##
## Generate ERA5 time series from a set of lat-longs using SSMT weather tools
##############################################################################################################

from simtools.AssetManager.FileList import FileList
from simtools.Managers.WorkItemManager import WorkItemManager
from simtools.SetupParser import SetupParser


def main():
    wi_name = 'ERA5 weather generation'

    # A .csv or a demographics file containing input coordinates
    root_dir = "/Users/bertozzivill/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/results/v3_era5_climate_rescaled/africa/02_kmeans"
    points_file = "site_ids_9_cluster.csv"

    # Start/end dates in one of the formats: year (2015) or year and day-of-year (2015032) or year-month-day (20150201)
    start_date = 2014
    end_date = 2015

    # Optional arguments


    # Data source selection: use "--ds" to select a data source. Currently available data sources are:
    #   ERA5: world-wide daily estimates for air/ground temperature, humidity and rainfall at 30km resolution.
    #   TAMSATv3: Daily rainfall estimates for all of Africa at 4km resolution.

    # CSV output format: use "--outtype csvfile" to generate weather files in .csv format.
    # Consider the size of the output CSV file. For example, for 1 node and 1 year the size will be ~20KB.

    optional_args = '--ds ERA5 --id-ref "Custom user" --node-col cluster --outtype csvfile'

    # To run a specific version add a tag (for example, "weather-files:1.1").
    # See available versions here: https://github.com/InstituteforDiseaseModeling/dst-era5-weather-data-tools/releases.
    docker_image = "weather-files"
    command_pattern = "python /app/generate_weather_asset_collection.py {} {} {} {}"
    command = command_pattern.format(points_file, start_date, end_date, optional_args)
    user_files = FileList(root=root_dir, files_in_root=[points_file])

    wi = WorkItemManager(item_name=wi_name, docker_image=docker_image, command=command, user_files=user_files,
        tags={'Demo': 'dtk-tools Docker WorkItem', 'WorkItem type': 'Docker', 'Command': command })

    wi.execute()


if __name__ == "__main__":
    SetupParser.default_block = 'HPC'
    SetupParser.init()
    main()


