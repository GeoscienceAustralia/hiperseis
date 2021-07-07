from pprint import pprint
import requests

def get_list_of_events(start_time, end_time):
    url = "https://earthquakes.ga.gov.au/geoserver/earthquakes/wfs?service=WFS&request=getfeature&typeNames=earthquakes:earthquakes&outputFormat=application/json&CQL_FILTER=display_flag=%27Y%27%20AND%20located_in_australia=%27Y%27%20AND%20preferred_magnitude%3E=0%20AND%20preferred_magnitude%3C=9.94%20AND%20origin_time%20BETWEEN%20"+start_time+"Z%20AND%20"+end_time+"Z"
    response = requests.get(url)
    return response.json()

def get_earthquakes_focal_mechanism(earthquake_id):
    url = "https://earthquakes.ga.gov.au/geoserver/earthquakes/wfs?service=WFS&request=getfeature&typeNames=earthquakes:earthquakes_focal_mechanism&outputFormat=application/json&CQL_FILTER=earthquake_id="+str(earthquake_id)
    response = requests.get(url)
    return response.json()

def get_earthquakes_shakemap(event_id):
    url = "https://earthquakes.ga.gov.au/geoserver/earthquakes/wfs?service=WFS&request=getfeature&typeNames=earthquakes:shakemap&outputFormat=application/json&CQL_FILTER=event_id=%27"+str(event_id)+"%27&dt=1625543755701"
    response = requests.get(url)
    return response.json()

def get_felt_reports(event_id):
    url = "https://earthquakes.ga.gov.au/geoserver/earthquakes/wfs?service=WFS&request=getfeature&typeNames=earthquakes:earthquakes_felt_reports&outputFormat=application/json&CQL_FILTER=event_id="+"'"+event_id+"'"
    response = requests.get(url)
    return response.json()

def get_trace_diagram(event_id):
    url = "https://cdn.eatws.net/skip/events/"+event_id+"/traces.png"
    open('traces.png', 'wb').write(requests.get(url, allow_redirects=True).content)
    print('file saved: traces.png')

def get_station_information(earthquake_id):
    url = "https://earthquakes.ga.gov.au/geoserver/earthquakes/wfs?service=WFS&request=getfeature&typeNames=earthquakes:earthquakes_stations&outputFormat=application/json&CQL_FILTER=earthquake_id="+str(earthquake_id)
    response = requests.get(url)
    return response.json()

def get_magnitutes_information(earthquake_id):
    url = "https://earthquakes.ga.gov.au/geoserver/earthquakes/wfs?service=WFS&request=getfeature&typeNames=earthquakes:earthquakes_magnitude&outputFormat=application/json&CQL_FILTER=earthquake_id="+str(earthquake_id)
    response = requests.get(url)
    return response.json()

def main():
    start_time = "2020-06-24T09:27:00"
    end_time = "2021-06-24T09:27:00"
    list_of_events = get_list_of_events(start_time, end_time)["features"]

    for event_details in list_of_events:

        print("event_details")
        pprint(event_details)

        earthquake_id = event_details["properties"]["earthquake_id"]
        event_id = event_details["properties"]["event_id"]

        focal_mechanism_information = get_earthquakes_focal_mechanism(earthquake_id)
        print("focal_mechanism_information")
        pprint(focal_mechanism_information)


        shakemap_information = get_earthquakes_shakemap(event_id)
        print("shakemap_information")
        pprint(shakemap_information)

        felt_reports = get_felt_reports(event_id)
        print("felt_reports")
        pprint(felt_reports)

        print("trace_diagram")
        get_trace_diagram(event_id)

        print("station_information")
        station_information = get_station_information(earthquake_id)
        pprint(station_information)

        print("magnitutes_information")
        magnitutes_information = get_magnitutes_information(earthquake_id)
        pprint(magnitutes_information)
        break

if __name__ == '__main__':
    main()
