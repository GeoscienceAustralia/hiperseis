# coding: utf-8

# # Visualization of events, station, rays in csv files 
# 
# 
# ##  Fei Zhang
# ### 2018-09-06
#
# 
# -  1. Python Pandas plot function (current cluster script)
# -  2. GMT python interface
# -  3. Geopandas mapping 
# 
# 
# ## installation of gmt-6 and python interface gmt-python
# 
# - conda install gmt -c conda-forge/label/dev
# - conda install numpy pandas xarray
# - pip install https://github.com/GenericMappingTools/gmt-python/archive/master.zip
# 
# 
# installation of geopandas
#   conda install geopandas

import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
# import sys
from shapely.geometry import mapping
from shapely.geometry import Point
from shapely.geometry import Point, Polygon, LineString, LinearRing

# sys.path.append("/Softlab/Githubz/passive-seismic")
# from seismic.cluster.cluster import Grid2
from seismic.traveltime.sort_rays import sort, sort2


def read_csv2pdf(csvfile):
    """
    Read in a csv file into a pandas dataframe. 
    Make sure the column names match the csv files. 
    delimiter/separator is whitespace or comma
    """

    # read infile, tweek below
    in_pdf = pd.read_csv(csvfile)  # assume there is header line=0; separator = comma,
    # finalpdf =  pd.read_csv(csvfile, header=None,  names=col_names ) #  no header line, separator = comma,
    # finalpdf =  pd.read_csv(csvfile,  sep='\s+', header=None,  names=col_names) # white space

    # columns you are interested?
    col_names = ['source_block', 'station_block', 'residual', 'event_number',
                 'source_longitude', 'source_latitude', 'source_depth',
                 'station_longitude', 'station_latitude', 'observed_tt', 'locations2degrees', 'station_code', 'SNR',
                 'P_or_S']

    # finalpdf = in_pdf[col_names]

    return in_pdf  # pandas_data_frame


def sort_plus(inputcsv=None):
    # P wave
    # inputcsv='/g/data/ha3/fxz547/Githubz/passive-seismic/seismic_events_arrivals_P_0.csv'
    # inputcsv='/g/data/ha3/fxz547/Githubz/passive-seismic/tempworks/outfile_P.csv'
    inputcsv = '/g/data/ha3/fxz547/travel_time_tomography/CSV_NewFormatAug10/FZ01-pst-cluster2/run3/P_out.csv'
    inputcsv = '/g/data/ha3/fxz547/travel_time_tomography/CSV_NewFormatAug10/FZ01-pst-cluster2/run3/P_out_translated.csv'
    # inputcsv='/Softlab/travel_time_tomography/PST/CSV_New_FZ01-pst-cluster2_run3/P_out.csv'
    residual_cutoff = 5.0  # cutoff value for P is 5s
    sortedfile = 'sortedfile_P.csv'
    sortedfile2 = 'sortedfile2_P.csv'

    # S wave
    # inputcsv = '/g/data1a/ha3/fxz547/travel_time_tomography/run5_events_1deg/outfile_S.csv'
    # inputcsv = '/g/data/ha3/fxz547/travel_time_tomography/CSV_NewFormatAug10/FZ01-pst-cluster2/run3/S_out.csv'
    # residual_cutoff = 10.0
    # sortedfile = 'sortedfile_S.csv'
    # sortedfile2 = 'sortedfile2_S.csv'

    mypdf = sort(inputcsv, sortedfile, residual_cutoff)  # select the median travel time
    # finalpdf2 = sort2(inputcsv,sortedfile,residual_cutoff)  # select the best SNR ray

    # show property of mypdf

    print(mypdf.shape)

    print(mypdf.head())

    print(mypdf.source_block.nunique())  # number of unique values

    # number of unique station_block number.
    # This may be less than the unique station_code, due to some nearby stations are clustered into one.
    print(mypdf.station_block.nunique())  # number of unique values.

    # number of unique stations
    print(mypdf.station_code.nunique())

    print(mypdf.groupby(['station_block', 'station_code']).count().shape)

    return mypdf  # sorted reduced


def plot1(mypdf):
    """
    :param mypdf: which dataframe to view in the following??
    :return:
    """

    # plt.figure(); mypdf.plot(x='event_number', y='observed_tt')
    # plt.figure(); mypdf.plot(x='event_number', y='source_depth')
    # plt.figure(); mypdf.plot(x='event_number', y='locations2degrees')

    # plt.figure(); mypdf.plot.scatter(x='event_number', y='residual',figsize=(12,8))  # less than +-10s

    # plt.figure(); mypdf.plot.scatter(x='event_number', y='source_block',figsize=(12,8))

    # plt.figure(); mypdf.plot.scatter(x='source_block', y='station_block',figsize=(12,8))

    # # 1. Python Pandas package plot function (with basemap)

    mypdf.plot.scatter(x='source_longitude', y='source_latitude', figsize=(12, 8))
    plt.show()

    mypdf.plot.scatter(x='station_longitude', y='station_latitude', figsize=(12, 8))
    plt.show()

    mypdf.plot.scatter(x='locations2degrees', y='observed_tt', figsize=(12, 8))

    plt.show()

    return


def plot_gmt(mypdf):
    """
    GMT-python plotting map figure saved into PNG file.
    But not displaying on screen

    Require to conda install gmt-python module in anaconda python
    gmt version-6 will be installed.
    :param mypdf:
    :return:
    """

    import gmt

    my_region = [mypdf.source_longitude.min(), mypdf.source_longitude.max(),
                 mypdf.source_latitude.min() - 1, mypdf.source_latitude.max() + 1]

    print(my_region)

    fig = gmt.Figure()
    fig.coast(region=my_region, projection='M6i', frame=True,
              land='black', water='skyblue')
    fig.plot(x=mypdf.source_longitude, y=mypdf.source_latitude,
             style='c0.3c', color='white', pen='black')

    fig.savefig("events.png")

    fig.show(dpi=400, width=1000)
    plt.show()

    fig = gmt.Figure()
    fig.coast(region=my_region, projection='M6i', frame=True, land='black', water='skyblue')
    fig.plot(x=mypdf.station_longitude, y=mypdf.station_latitude,
             style='c0.3c', color='white', pen='black')

    fig.savefig("station.png")
    fig.show(dpi=400, width=1000)
    plt.show()


def plot_geopandas(mypdf):
    """
    Geopandas and Maps
    :param mypdf:
    :return:
    """


    event_locations = [Point(xy) for xy in zip(mypdf.source_longitude, mypdf.source_latitude)]
    # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
    # if you want to df = df.drop(['Lon', 'Lat'], axis=1)
    mycrs = {'init': 'epsg:4326'}  # WGS84
    mycrs = {'init': 'epsg:4283'}  # GDA94
    geopdf = gpd.GeoDataFrame(mypdf, crs=mycrs, geometry=event_locations)

    # myax = geopdf.plot(figsize=[20,10])

    # myax.set_xlabel('Longitude')
    # myax.set_ylabel('Latitude+
    # title_str= "event locations"
    # myax.set_title(title_str)

    # geopandas included shape datasets
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    mymap = world.plot(alpha=0.5, figsize=(20, 10))

    mymap.set_xlim([50, 180])
    mymap.set_ylim([-70, 70])

    mymap.set_xlabel('Longitude')
    mymap.set_ylabel('Latitude')
    title_str = "Seismic events locations on Map"
    mymap.set_title(title_str)
    geopdf.plot(ax=mymap, marker='o', color='red', markersize=2)

    plt.savefig('event_locations_on_map.png')

    plt.show()
    return


def plot_rays(mypdf):
    """

    :param mypdf:
    :return:
    """
    mycrs = {'init': 'epsg:4326'}  # WGS84

    # geo_df = gpd.GeoDataFrame(pdf, crs=crs, geometry=mt_locations)

    mypdf['ray'] = mypdf.apply(lambda x: LineString([(x.source_longitude, x.source_latitude),
                                                     (x.station_longitude, x.station_latitude)]), axis=1)

    geopdf_ray = gpd.GeoDataFrame(mypdf, crs=mycrs, geometry='ray')

    # myax=geopdf_ray.plot(figsize=[20,10])

    # myax.set_xlabel('Longitude')
    # myax.set_ylabel('Latitude')
    # title_str= "event->station Rays"
    # myax.set_title(title_str)

    # geopandas included shape datasets
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    mymap = world.plot(alpha=0.5, figsize=(20, 10))

    # mymap.set_xlim([-180,180])
    # mymap.set_ylim([-80,80])

    mymap.set_xlim([50, 180])
    mymap.set_ylim([-70, 70])

    mymap.set_xlabel('Longitude')
    mymap.set_ylabel('Latitude')
    title_str = "event->station Rays"
    mymap.set_title(title_str)

    geopdf_ray.plot(ax=mymap)
    plt.show()

    return


if __name__ == "__main__":
    newpdf = sort_plus()

    #plot1(newpdf)

    plot_geopandas(newpdf)

    #plot_gmt(newpdf)
