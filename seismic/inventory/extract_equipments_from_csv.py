#!/usr/bin/env python
# coding: utf-8
"""
Extract Sensor and Digitizer Equipment Information from a CSV file

Fei Zhang
2020 Sept.

Input file is provided /g/data1a/ha3/fxz547/Githubz/hiperseis/tests/testdata/FieldSiteVisitLive.csv (4781 rows, 36 columns)
Algorithm: use pandas to read and subset the csv file. filter the rows, select the required info, output to a new csv file
"""

import pandas as pd
import json


class EquipmentExtractor:

    def __init__(self, csvfile="../../tests/testdata/FieldSiteVisitLive.csv"):
        self.csvfile = csvfile
        # where is the input CSV file

        self.equipments = self.get_equipment_from_csv()

    def get_equipment_from_csv(self):
        pdf = pd.read_csv(self.csvfile, delimiter=',')

        my_columns = ["submissionid",
                      "submissiondatetime",
                      "deviceid",
                      "username",
                      "Select_Site_ID",
                      "GPS_Locations_Subform.GPS_Location_Group.Select_Location_Type",
                      "GPS_Locations_Subform.GPS_Location_Group.Set_Location",
                      "Replace_Equipment_Yes_No",
                      "Asset_Tracking_Subform.Asset_Barcodes_Group.DB_Asset_Type",
                      "Asset_Tracking_Subform.Asset_Barcodes_Group.DB_Asset_UniqueID"]

        mypdf = pdf[my_columns]

        my_criteria = ((mypdf[
                            'Asset_Tracking_Subform.Asset_Barcodes_Group.DB_Asset_Type'] == "Nanometrics Trillium Compact 120s") |
                       (mypdf['Asset_Tracking_Subform.Asset_Barcodes_Group.DB_Asset_Type'] == "Guralp Minimus"))

        # my_criteria

        my_selected_pdf = mypdf.loc[my_criteria].copy()  # use .copy() to fix "SettingWithCopyWarning"

        my_selected_pdf.head()

        my_selected_pdf.shape

        output_cols = ["Select_Site_ID",
                       "Asset_Tracking_Subform.Asset_Barcodes_Group.DB_Asset_Type",
                       "Asset_Tracking_Subform.Asset_Barcodes_Group.DB_Asset_UniqueID"]

        new_col_dict ={
            "Select_Site_ID": "NetSta",
            "Asset_Tracking_Subform.Asset_Barcodes_Group.DB_Asset_Type":"Description",
            "Asset_Tracking_Subform.Asset_Barcodes_Group.DB_Asset_UniqueID":"SerNumber"
        }

        return my_selected_pdf[output_cols].rename(columns = new_col_dict)

    def get_sensors(self, net, sta):
        """ get the sensors for a specific (net,sta)-code pair
        return a json object
        """

        _criteria = (self.equipments["Description"] == "Nanometrics Trillium Compact 120s") & \
                    (self.equipments["NetSta"] == str(net + str(sta)))

        _equip = self.equipments.loc[_criteria].copy()

        json_str = _equip.to_json(orient="records")
        return json.loads(json_str)

    def get_digitizer(self, net, sta):
        """ get the digitizers for a specific (net,sta)-code pair
         return a json object
        """

        _criteria = (self.equipments["Description"] == "Guralp Minimus") & \
                    (self.equipments["NetSta"] == str(net + str(sta)))

        _equip=self.equipments.loc[_criteria].copy()

        json_str = _equip.to_json(orient="records")
        return json.loads(json_str)

# =======================
# A quick test run
if __name__ == '__main__':
    my_equip_obj = EquipmentExtractor() #provide the original csv files from field work company

    pdf_equip = my_equip_obj.get_equipment_from_csv()

    pdf_equip.to_csv("OA_Equipments.csv", index=False)

    # OA station codes
    station_codes=["BS24", "BS25", "BS26","BS27",  "BS28", "BT23", "BT24", "BT25", "BT26","BT27", "BT28"]

    for stacode in station_codes:
        # print("Digitizer: ", type(my_equip_obj.get_digitizer("OA", stacode)))  # list

        if len(my_equip_obj.get_sensors("OA", stacode))>0:
            print("Sensor: ",  my_equip_obj.get_sensors("OA", stacode)[0].get('Description'))
        else:
            print("%s %s No sensors" % ("OA", stacode))


