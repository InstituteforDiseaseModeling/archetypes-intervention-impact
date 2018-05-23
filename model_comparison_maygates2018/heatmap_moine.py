
# heatmap 1: just pixels under hfs
print("extracting hf covariate")

cov_path = os.path.join(out_path, "rasters", input_name)
files = os.listdir(cov_path)
all_area_tifs = list(filter(re.compile(r'{name}_all.*tif'.format(name=input_name)).search, files))
all_area_tifs.sort()

all_shp = make_shapefile(hf_data.copy(), type='point',
                         to_crs=shp.crs, lat_name='lat', lon_name='lng')
shapes = [shapely.geometry.mapping(g) for g in all_shp['geometry']]

for idx, fname in enumerate(all_area_tifs):
    month = int(re.match('.*\.([0-9]{2})\.', fname).groups()[0])
    full_fname = os.path.join(cov_path, fname)
    all_shp['{name}_{month}'.format(name=input_name, month=month)] = extract_latlongs(full_fname, shapes)

all_shp.sort_values(by=['hf'], inplace=True)
all_shp = all_shp.query('hf!="Caputine"')
to_plot = all_shp.filter(regex=input_name, axis=1)
to_plot.columns = list(range(1, 13))
to_plot.columns.name = 'month'
to_plot.index = all_shp.hf
to_plot.index.name = "hf_name"

fig, axes = plt.subplots(1, 2, figsize=(14, 7))
sns.heatmap(to_plot, cmap="BuPu", ax=axes[0])
axes[0].set_title("{name}, Magude HFs".format(name=input_name))
axes[0].set_ylabel("")

to_plot.to_csv(os.path.join(cov_path, "hf_{name}.csv".format(name=input_name)))


# heatmap 2: mean values of all pixels within catchment
print("finding means over full catchment areas")


def find_mean(name):
    # todo: modify this to actually remove the appropriate "nodata" value
    in_raster = rasterio.open(name)
    in_data = in_raster.read(1)
    mean = in_data[in_data >= 0].mean()

    return mean


mean_hfca_cov = pd.DataFrame({'path': list(filter(re.compile(r'.*tif').search, files))})
mean_hfca_cov['month'] = mean_hfca_cov['path'].apply(lambda x: int(re.match('.*\.([0-9]{2})\.', x).groups()[0]))
mean_hfca_cov['hf_name'] = mean_hfca_cov['path'].apply(lambda x: re.match('{name}_(.+)\.Month'.format(name=input_name), x).groups()[0])
mean_hfca_cov['full_path'] = mean_hfca_cov['path'].apply(lambda x: os.path.join(cov_path, x))
mean_hfca_cov['mean_cov'] = mean_hfca_cov['full_path'].apply(find_mean)

to_plot = pd.pivot_table(mean_hfca_cov, values='mean_cov', columns='month', index=['hf_name'])

sns.heatmap(to_plot.drop('all'), cmap="BuPu", ax=axes[1])
axes[1].set_title("{name}, Mean Magude HFCAs".format(name=input_name))
axes[1].set_ylabel('')
plt.savefig(os.path.join(cov_path, "{name}_heat.jpg".format(name=input_name)))

to_plot.to_csv(os.path.join(cov_path, "hfca_mean_{name}.csv".format(name=input_name)))

# line plot of mean time series
line_plot = to_plot.reset_index().melt(id_vars="hf_name", value_name=input_name)

fig, axes = plt.subplots(1, 2, figsize=(14, 7))
sns.pointplot(x="month", y=input_name, hue="hf_name", data=line_plot, ax=axes[0])
sns.pointplot(x="month", y=input_name, data=line_plot.query('hf_name=="Moine"'), ax=axes[1])
axes[1].set_title("Mean {name}, Moine".format(name=input_name))

plt.savefig(os.path.join(cov_path, "mean_{name}_line.jpg".format(name=input_name)))

# cluster by time
print("clustering mean cov")
time = line_plot.query('hf_name=="Moine"')[['month', input_name]]

for i, n in enumerate([2, 3, 4]):

    for_clust = time[input_name].as_matrix().reshape(-1, 1)
    clust = KMeans(n_clusters=n).fit(for_clust)
    col_name = 'clust_{count}'.format(count=n)
    time[col_name] = clust.labels_

time = time.melt(id_vars=['month', input_name])
g = sns.FacetGrid(time, col='variable', hue='value')
g = (g.map(plt.scatter, "month", input_name))

plt.savefig(os.path.join(cov_path, "cluster_Moine_{name}.jpg".format(name=input_name)))


# focus on Moine: heatmaps under households
print("Household Level Heatmap, Moine HFC")
moine_shp = make_shapefile(hh_data.query('hf_name=="Moine"').copy(), type='point',
                           to_crs=shp.crs, lat_name='lat_r2', lon_name='lng_r2')

shapes = [shapely.geometry.mapping(g) for g in moine_shp['geometry']]

for idx, fname in enumerate(all_area_tifs):
    month = int(re.match('.*\.([0-9]{2})\.', fname).groups()[0])
    full_fname = os.path.join(cov_path, fname)
    moine_shp['{name}_{month}'.format(name=input_name, month=month)] = extract_latlongs(full_fname, shapes)

moine_shp['latitude'] = moine_shp['geometry'].apply(lambda geom: geom.x)
moine_shp.sort_values(by=['latitude'], inplace=True)
to_plot = moine_shp.filter(regex=input_name, axis=1)
to_plot.columns = list(range(1, 13))
to_plot.columns.name = 'month'

fig, axes = plt.subplots(1,2, figsize=(14, 7))
sns.heatmap(to_plot, cmap="viridis", ax=axes[0])
axes[0].set_title("{name}, Moine Households".format(name=input_name))
axes[0].set_ylabel("")

# and heatmap under every pixel
moine_tifs = list(filter(re.compile(r'{name}_Moine.*tif'.format(name=input_name)).search, files))
moine_tifs.sort()

pixels = {}

for idx, fname in enumerate(moine_tifs):
    month = int(re.match('.*\.([0-9]{2})\.', fname).groups()[0])
    full_fname = os.path.join(cov_path, fname)

    raster_data = rasterio.open(full_fname).read(1)

    pixels[month] = raster_data[raster_data >= 0]

pixels = pd.DataFrame(pixels)
pixels.columns.name = "month"
sns.heatmap(pixels, cmap="viridis", ax=axes[1])
axes[1].set_title("{name}, Moine Pixels".format(name=input_name))
axes[1].set_ylabel("")

plt.savefig(os.path.join(cov_path, "moine_heat.jpg"))
