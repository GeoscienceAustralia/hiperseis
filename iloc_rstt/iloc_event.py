from obspy.core.event import Event, Catalog, Origin, Comment
from obspy import read_events
from inventory.parse_inventory import read_all_stations
from iloc_rstt.config import event_file


class ILocEvent:
    """
    Just a convenience class for event enhancement using iLoc
    """
    def __init__(self, event):
        self.event = event

    def update(self):
        self.event.picks += self.add_picks()
        self.event.origins += self.add_origin()

    def add_picks(self):
        pass

    def add_origin(self):
        return Origin()

    def add_stations(self):
        self._add_primary_stations()
        self._add_temporary_stations()

    def _add_primary_stations(self):
        pass

    def _add_temporary_stations(self):
        pass


class ILocCatalog:
    """
    Class that supports event/catalog enhancement using iLoc.
    iLoc ref: http://www.seismology.hu/index.php/en/home/iloc
    """
    def __init__(self, event_xml):
        self.event_xml = event_xml
        self.cat = read_events(event_xml)
        self.events = self.cat.events
        self.iloc_events = []

    def update(self):
        for e in self.events:
            self.iloc_events.append(ILocEvent(e).update())

    def write(self, new_sc3ml):
        self.cat.write(new_sc3ml)


if __name__ == "__main__":
    cat = ILocCatalog(event_file)
    cat.update()
