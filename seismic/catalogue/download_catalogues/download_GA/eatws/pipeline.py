from pathlib import Path
from .acquire_data import acquire_native_eatws_data
from .format_data import save_eatws_data_to_quakeml


def download_eatws_event_data(start_time, end_time, min_magnitude, max_magnitude, bounding_box, output_directory):
    events_data = acquire_native_eatws_data(start_time, end_time, min_magnitude, max_magnitude, bounding_box)
    print("Total events found: " + str(len(events_data)))
    for i, event_data in enumerate(events_data):
        print("Currently processing: " + str(i) + "/" + str(len(events_data)))
        origin_time = event_data["event_details"]["properties"]["origin_time"]
        print("Event origin time: " + str(origin_time))
        origin_time = origin_time.replace(":", "-")
        output_event_file = Path(output_directory) / Path(origin_time + ".xml")
        save_eatws_data_to_quakeml(event_data, output_event_file)
