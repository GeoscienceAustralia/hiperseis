import json
from shapely.geometry.polygon import Polygon, Point

class StandardDomain():
    def __init__(self, domain_json_file:str):
        dom = json.load(open(domain_json_file, 'r'))
        lat_center, lat_extent, lon_center, lon_extent = \
            dom['lat_lon_box']['latitude_center'], \
            dom['lat_lon_box']['latitude_extent'], \
            dom['lat_lon_box']['longitude_center'], \
            dom['lat_lon_box']['longitude_extent']

        c = [lon_center, lat_center]
        e = [lon_extent, lat_extent]

        coords = ((c[0] - e[0], c[1] + e[1]),
                  (c[0] - e[0], c[1] - e[1]),
                  ((c[0] + e[0]), c[1] - e[1]),
                  ((c[0] + e[0]), c[1] + e[1]))
        self.coords = np.array(coords)
        self.bounding_polygon = Polygon(coords)
    # end func

    def contains(self, lon:float, lat:float):
        p = Point((lon, lat))

        return self.bounding_polygon.contains(p)
    # end func
# end class

class Domain():
    def __init__(self, domain_json_file:str):
        dom = json.load(open(domain_json_file, 'r'))

        coords = dom['features'][0]['geometry']['coordinates'][0]
        self.bounding_polygon = Polygon(coords)
    # end func

    def contains(self, lon:float, lat:float):
        if(lon < 0): lon += 360 # temporary hack to cater for dateline crossings
        p = Point((lon, lat))

        return self.bounding_polygon.contains(p)
    # end func
# end class
