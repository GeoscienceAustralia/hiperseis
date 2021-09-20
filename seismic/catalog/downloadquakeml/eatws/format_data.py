import pickle
from obspy import UTCDateTime
from obspy.core.event import Magnitude, Catalog, Event, Origin, Pick, Arrival, CreationInfo, WaveformStreamID, \
    QuantityError, OriginQuality, OriginUncertainty, EventDescription
from obspy.geodetics import FlinnEngdahl
from pprint import pprint


def calculate_backazimuth_from_azimuth(azimuth):
    if azimuth < 180:
        return azimuth + 180
    return azimuth - 180


def get_arrivals_and_picks(event_data):
    arrivals = []
    picks = []
    for station_information in event_data["station_information"]["features"]:
        pick_object = Pick(
            resource_id="pick_id_" + str(station_information["properties"]["arrival_id"]),
            force_resource_id=True,
            time=UTCDateTime(station_information["properties"]["arrival_time"]),
            time_errors=QuantityError(),
            waveform_id=WaveformStreamID(
                network_code=station_information["properties"]["network_code"],
                station_code=station_information["properties"]["station_code"],
                location_code="",
                channel_code=station_information["properties"]["channel_code"]
            ),
            # filter_id=None,
            # method_id=None,
            # horizontal_slowness=None,
            # horizontal_slowness_errors=None,
            backazimuth=calculate_backazimuth_from_azimuth(station_information["properties"]["azimuth"]),
            backazimuth_errors=QuantityError(),
            # slowness_method_id=None,
            # onset=None,
            phase_hint=station_information["properties"]["phase"],
            polarity="undecidable",
            evaluation_mode=event_data["event_details"]["properties"]["evaluation_mode"],
            evaluation_status=event_data["event_details"]["properties"]["evaluation_status"],
            comments=[],
            creation_info=CreationInfo(
                agency_id=event_data["event_details"]["properties"]["source"],
                author=event_data["event_details"]["properties"]["source"],
                author_uri="ga.gov.au",
                version="1.0"
            )
        )

        arrival_object = Arrival(
            resource_id="arrival_id_" + str(station_information["properties"]["arrival_id"]),
            force_resource_id=False,
            pick_id="pick_id_" + str(station_information["properties"]["arrival_id"]),
            phase=station_information["properties"]["phase"],
            # time_correction=None,
            azimuth=station_information["properties"]["azimuth"],
            distance=station_information["properties"]["distance"],
            # takeoff_angle=None,
            takeoff_angle_errors=QuantityError(),
            # time_residual=None,
            # horizontal_slowness_residual=None,
            # backazimuth_residual=None,
            # time_weight=None,
            # horizontal_slowness_weight=None,
            # backazimuth_weight=None,
            earth_model_id=event_data["event_details"]["properties"]["earth_model_id"],
            # comments=[],
            creation_info=CreationInfo(
                agency_id=event_data["event_details"]["properties"]["source"],
                author=event_data["event_details"]["properties"]["source"],
                author_uri="ga.gov.au",
                version="1.0"
            )
        )
        arrivals.append(arrival_object)
        picks.append(pick_object)

    return arrivals, picks


def get_origins(event_data, arrivals, origin_id):
    origins = [Origin(
        resource_id=origin_id,
        force_resource_id=False,
        time=UTCDateTime(event_data["event_details"]["properties"]["origin_time"]),
        time_errors=QuantityError(
            uncertainty=event_data["event_details"]["properties"]["origin_time_uncertainty"]
        ),
        longitude=event_data["event_details"]["properties"]["latitude"],
        longitude_errors=QuantityError(),
        latitude=event_data["event_details"]["properties"]["longitude"],
        latitude_errors=QuantityError(),
        depth=event_data["event_details"]["properties"]["depth"],
        depth_errors=QuantityError(
            uncertainty=event_data["event_details"]["properties"]["depth_uncertainty"]
        ),
        depth_type="other",
        # time_fixed=None,
        # epicenter_fixed=None,
        # reference_system_id=None,
        # method_id=None,
        earth_model_id=event_data["event_details"]["properties"]["earth_model_id"],
        arrivals=arrivals,
        # composite_times=None,
        quality=OriginQuality(
            # associated_phase_count=None,
            used_phase_count=event_data["event_details"]["properties"]["phase_count"],
            # associated_station_count=None,
            used_station_count=event_data["event_details"]["properties"]["station_count"],
            # depth_phase_count=None,
            standard_error=event_data["event_details"]["properties"]["standard_error"],
            azimuthal_gap=event_data["event_details"]["properties"]["azimuthal_gap"],
            # secondary_azimuthal_gap=None,
            # ground_truth_level=None,
            minimum_distance=event_data["event_details"]["properties"]["minimum_distance"],
            maximum_distance=event_data["event_details"]["properties"]["maximum_distance"],
            # median_distance=None,
        ),
        # origin_type=None,
        origin_uncertainty=OriginUncertainty(
            # horizontal_uncertainty=None,
            min_horizontal_uncertainty=event_data["event_details"]["properties"]["min_horizontal_uncertainty"],
            max_horizontal_uncertainty=event_data["event_details"]["properties"]["max_horizontal_uncertainty"],
            azimuth_max_horizontal_uncertainty=event_data["event_details"]["properties"][
                "azimuth_horizontal_uncertainty"],
            # confidence_ellipsoid=None,
            # preferred_description=None,
            # confidence_level=None
        ),
        region=event_data["station_information"]["features"][0]["properties"]["network_code"],
        # evaluation_mode=None,
        # evaluation_status=None,
        # comments=[],
        creation_info=CreationInfo(
            agency_id=event_data["event_details"]["properties"]["source"],
            author=event_data["event_details"]["properties"]["source"],
            author_uri="ga.gov.au",
            version="1.0"
        )
    )]
    return origins


def get_magnitude(event_data, origin_id):
    for magnitude_information in event_data["magnitudes_information"]["features"]:
        if magnitude_information["properties"]["evaluation_status"] == "confirmed":
            magnitudes = [Magnitude(
                resource_id="magnitude_id_" + str(magnitude_information["properties"]["earthquake_id"]),
                force_resource_id=False,
                mag=magnitude_information["properties"]["magnitude"],
                mag_errors=QuantityError(),
                magnitude_type=magnitude_information["properties"]["type"],
                origin_id=origin_id,
                # method_id=None,
                station_count=event_data["event_details"]["properties"]["station_count"],
                azimuthal_gap=event_data["event_details"]["properties"]["azimuthal_gap"],
                evaluation_mode=event_data["event_details"]["properties"]["evaluation_mode"],
                evaluation_status=event_data["event_details"]["properties"]["evaluation_status"],
                comments=[],
                station_magnitude_contributions="",
                creation_info=CreationInfo(
                    agency_id=event_data["event_details"]["properties"]["source"],
                    author=event_data["event_details"]["properties"]["source"],
                    author_uri="ga.gov.au",
                    version="1.0"
                )
            )]
            return magnitudes


def get_event(event_data, picks, origins, magnitudes):
    event = Event(
        resource_id="resource_id_" + event_data["event_details"]["properties"]["event_id"],
        force_resource_id=False,
        event_type="earthquake",
        event_type_certainty="known",
        creation_info=CreationInfo(
            agency_id=event_data["event_details"]["properties"]["source"],
            author=event_data["event_details"]["properties"]["source"],
            author_uri="ga.gov.au",
            version="1.0"
        ),
        # event_descriptions=event_data["event_details"]["properties"]["description"],
        comments=[],
        picks=picks,
        # amplitudes=None,
        # focal_mechanisms=None,
        origins=origins,
        magnitudes=magnitudes,
        # station_magnitudes=None
    )
    return event


def save_eatws_data_to_quakeml(event_data, output_data_file):
    arrivals, picks = get_arrivals_and_picks(event_data)
    origin_id = "origins_id_" + event_data["event_details"]["properties"]["origin_id"]
    origins = get_origins(event_data, arrivals, origin_id)
    magnitudes = get_magnitude(event_data, origin_id)
    event = get_event(event_data, picks, origins, magnitudes)

    catalog = Catalog(
        events=[event],
        resource_id="catalog_id_0",
        description="catalog of ga events",
        comments=[],
        creation_info=CreationInfo(
            agency_id=event_data["event_details"]["properties"]["source"],
            author=event_data["event_details"]["properties"]["source"],
            author_uri="ga.gov.au",
            version="1.0"
        )
    )
    catalog.write(output_data_file, format="QUAKEML")

    # Debugging
    # print(FlinnEngdahl().get_region(event_data["event_details"]["properties"]["longitude"], event_data["event_details"]["properties"]["latitude"]))


if __name__ == '__main__':
    input_data_file_local = 'C:/Users/sheec/Desktop/Project/genquakeml/data/test_event.pkl'
    output_data_file_local = "C:/Users/sheec/Desktop/Project/genquakeml/data/quakeml_events.xml"
    with open(input_data_file_local, 'rb') as handle:
        event_data_local = pickle.load(handle)

    save_eatws_data_to_quakeml(event_data_local, output_data_file_local)
