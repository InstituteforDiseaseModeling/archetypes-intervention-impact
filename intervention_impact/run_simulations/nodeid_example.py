
import pandas as pd
from dtk.tools.demographics.Node import Node, nodeid_from_lat_lon

res=30

sites = pd.read_csv("site_details.csv")
sites["facility_id"] = sites.apply(lambda row: nodeid_from_lat_lon(row["lat"], row["lon"], res), axis=1)

print(sites[["name", "lat", "lon", "facility_id"]])