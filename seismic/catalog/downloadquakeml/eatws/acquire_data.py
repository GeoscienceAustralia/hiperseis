import pickle
import requests
import logging
from pprint import pprint
from .constants import GA_EATWS_API


def get_list_of_events(start_time, end_time, min_magnitude, max_magnitude, bounding_box):
    url = GA_EATWS_API + "/wfs?service=WFS&request=getfeature&typeNames=earthquakes:earthquakes&outputFormat=application/json&CQL_FILTER=display_flag=%27Y%27%20AND%20located_in_australia=%27Y%27%20AND%20preferred_magnitude%3E=" + str(
        min_magnitude) + "%20AND%20preferred_magnitude%3C=" + str(
        max_magnitude) + "%20AND%20origin_time%20BETWEEN%20" + start_time + "Z%20AND%20" + end_time + "Z"
    print(url)
    response = requests.get(url)
    if response.status_code == 400:
        return {'features': [], 'totalFeatures': 0}

    # remove event if its not within the bounding_box
    min_latitude, min_longitude, max_latitude, max_longitude = bounding_box
    print("min_latitude, min_longitude, max_latitude, max_longitude", bounding_box)
    results = response.json()
    features = []
    for event in results["features"]:
        longitude, latitude = event["geometry"]["coordinates"]
        print("current event latitude, longitude", longitude, latitude)
        if min_longitude < longitude < max_longitude and min_latitude < latitude < max_latitude:
            features.append(event)
    results["features"] = features
    return results


def get_station_information(earthquake_id):
    url = GA_EATWS_API + "wfs?service=WFS&request=getfeature&typeNames=earthquakes:earthquakes_stations&outputFormat=application/json&CQL_FILTER=earthquake_id=" + str(
        earthquake_id)
    response = requests.get(url)
    return response.json()


def get_magnitudes_information(earthquake_id):
    url = GA_EATWS_API + "wfs?service=WFS&request=getfeature&typeNames=earthquakes:earthquakes_magnitude&outputFormat=application/json&CQL_FILTER=earthquake_id=" + str(
        earthquake_id)
    response = requests.get(url)
    return response.json()


def compile_event_information(event_details):
    event_data = {}

    event_data["event_details"] = event_details
    earthquake_id = event_details["properties"]["earthquake_id"]

    logging.info("Getting station information")
    station_information = get_station_information(earthquake_id)
    event_data["station_information"] = station_information

    logging.info("Getting magnitudes information")
    magnitudes_information = get_magnitudes_information(earthquake_id)
    event_data["magnitudes_information"] = magnitudes_information

    return event_data


def acquire_native_eatws_data(start_time, end_time, min_magnitude, max_magnitude, bounding_box):
    list_of_events = get_list_of_events(start_time, end_time, min_magnitude, max_magnitude, bounding_box)
    print("Total events found: " + str(list_of_events["totalFeatures"]))

    # Debugging
    # with open('output.txt', 'wt') as out:
    #     pprint(list_of_events, stream=out)

    events_data = []
    total_events = len(list_of_events["features"])
    for i, event_details in enumerate(list_of_events["features"]):
        print("Getting event details: " + str(i) + "/" + str(total_events))
        event_data = compile_event_information(event_details)
        events_data.append(event_data)
    return events_data


if __name__ == '__main__':
    event_data_local = acquire_native_eatws_data(start_time="2020-06-24T09:27:00", end_time="2021-06-24T09:27:00")

    with open('C:/Users/sheec/Desktop/Project/genquakeml/data/test_event.pkl', 'wb') as handle:
        pickle.dump(event_data_local, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Debugging
    # with open('data/test_event.pkl', 'rb') as handle:
    #     b = pickle.load(handle)
