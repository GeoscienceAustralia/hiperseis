#
#
# Gregion.pm: Flinn-Engdahl geographic region number functions
#
# Istvan Bondar, 2009/07/29
#
#
# Flinn-Engdahl regionalization hash
#     read Flinn-Engdahl regionalization file
#         $fe = read_FlinnEngdahl($fname)
# Geographic region number
#     Flinn-Engdahl geographic region number for a location
#         $grn = gregnum($lat, $lon, $fe);
# Geographic region name
#     Flinn-Engdahl geographic region name for a geographic region number
#         $grname = gregion($grn);
# Seismic region number
#     Flinn-Engdahl seismic region number for a geographic region number
#         $srn = gregtosreg($grn);
# Seismic region name
#     Flinn-Engdahl seismic region name for a seismic region number
#         $sname = sregion($srn);
#
# Flinn-Engdahl hash:
#     $fe->{error}  - set if error in reading file
#     $fe->{nlat}   - number of latitudes
#     @{$fe->{nl}}  - number of longitudes at each latitude(nlat)
#     @{$fe->{lon}} - array of longitudes(nlat x nl[i])
#     @{$fe->{grn}} - array of geographic region numbers(nlat x nl[i])
#
package Gregion;
require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(gregion gregtosreg sregion read_FlinnEngdahl gregnum gregn);
use strict;
use Interpolate;
use POSIX;

#
# $grname = gregion($grn)
#     returns Flinn-Engdahl geographic region name for a gregion number
#  Input Arguments:
#     $grn  - Flinn-Engdahl geographic region number
#  Return:
#     $grname - Flinn-Engdahl geographic region name
#
sub gregion {
    my $grn = shift;
    my @regname = (
        "Central Alaska",
        "Southern Alaska",
        "Bering Sea",
        "Komandorsky Islands region",
        "Near Islands",
        "Rat Islands",
        "Andreanof Islands",
        "Pribilof Islands",
        "Fox Islands",
        "Unimak Island region",
        "Bristol Bay",
        "Alaska Peninsula",
        "Kodiak Island region",
        "Kenai Peninsula",
        "Gulf of Alaska",
        "South of Aleutian Islands",
        "South of Alaska",
        "Southern Yukon Territory",
        "Southeastern Alaska",
        "Off coast of southeastern Alaska",
        "West of Vancouver Island",
        "Queen Charlotte Islands region",
        "British Columbia",
        "Alberta",
        "Vancouver Island region",
        "Off coast of Washington",
        "Near coast of Washington",
        "Washington-Oregon border region",
        "Washington",
        "Off coast of Oregon",
        "Near coast of Oregon",
        "Oregon",
        "Western Idaho",
        "Off coast of northern California",
        "Near coast of northern California",
        "Northern California",
        "Nevada",
        "Off coast of California",
        "Central California",
        "California-Nevada border region",
        "Southern Nevada",
        "Western Arizona",
        "Southern California",
        "California-Arizona border region",
        "California-Baja California border region",
        "Western Arizona-Sonora border region",
        "Off west coast of Baja California",
        "Baja California",
        "Gulf of California",
        "Sonora",
        "Off coast of central Mexico",
        "Near coast of central Mexico",
        "Revilla Gigedo Islands region",
        "Off coast of Jalisco",
        "Near coast of Jalisco",
        "Near coast of Michoacan",
        "Michoacan",
        "Near coast of Guerrero",
        "Guerrero",
        "Oaxaca",
        "Chiapas",
        "Mexico-Guatemala border region",
        "Off coast of Mexico",
        "Off coast of Michoacan",
        "Off coast of Guerrero",
        "Near coast of Oaxaca",
        "Off coast of Oaxaca",
        "Off coast of Chiapas",
        "Near coast of Chiapas",
        "Guatemala",
        "Near coast of Guatemala",
        "Honduras",
        "El Salvador",
        "Near coast of Nicaragua",
        "Nicaragua",
        "Off coast of central America",
        "Off coast of Costa Rica",
        "Costa Rica",
        "North of Panama",
        "Panama-Costa Rica border region",
        "Panama",
        "Panama-Colombia border region",
        "South of Panama",
        "Yucatan Peninsula",
        "Cuba region",
        "Jamaica region",
        "Haiti region",
        "Dominican Republic region",
        "Mona Passage",
        "Puerto Rico region",
        "Virgin Islands",
        "Leeward Islands",
        "Belize",
        "Caribbean Sea",
        "Windward Islands",
        "Near north coast of Colombia",
        "Near coast of Venezuela",
        "Trinidad",
        "Northern Colombia",
        "Lake Maracaibo",
        "Venezuela",
        "Near west coast of Colombia",
        "Colombia",
        "Off coast of Ecuador",
        "Near coast of Ecuador",
        "Colombia-Ecuador border region",
        "Ecuador",
        "Off coast of northern Peru",
        "Near coast of northern Peru",
        "Peru-Ecuador border region",
        "Northern Peru",
        "Peru-Brazil border region",
        "Western Brazil",
        "Off coast of Peru",
        "Near coast of Peru",
        "Central Peru",
        "Southern Peru",
        "Peru-Bolivia border region",
        "Northern Bolivia",
        "Central Bolivia",
        "Off coast of northern Chile",
        "Near coast of northern Chile",
        "Northern Chile",
        "Chile-Bolivia border region",
        "Southern Bolivia",
        "Paraguay",
        "Chile-Argentina border region",
        "Jujuy Province",
        "Salta Province",
        "Catamarca Province",
        "Tucuman Province",
        "Santiago del Estero Province",
        "Northeastern Argentina",
        "Off coast of central Chile",
        "Near coast of central Chile",
        "Central Chile",
        "San Juan Province",
        "La Rioja Province",
        "Mendoza Province",
        "San Luis Province",
        "Cordoba Province",
        "Uruguay",
        "Off coast of southern Chile",
        "Southern Chile",
        "Southern Chile-Argentina border region",
        "Southern Argentina",
        "Tierra del Fuego",
        "Falkland Islands region",
        "Drake Passage",
        "Scotia Sea",
        "South Georgia Island region",
        "South Georgia Rise",
        "South Sandwich Islands region",
        "South Shetland Islands",
        "Antarctic Peninsula",
        "Southwestern Atlantic Ocean",
        "Weddell Sea",
        "Off west coast of North Island",
        "North Island",
        "Off east coast of North Island",
        "Off west coast of South Island",
        "South Island",
        "Cook Strait",
        "Off east coast of South Island",
        "North of Macquarie Island",
        "Auckland Islands region",
        "Macquarie Island region",
        "South of New Zealand",
        "Samoa Islands region",
        "Samoa Islands",
        "South of Fiji Islands",
        "West of Tonga Islands (REGION NOT IN USE)",
        "Tonga Islands",
        "Tonga Islands region",
        "South of Tonga Islands",
        "North of New Zealand",
        "Kermadec Islands region",
        "Kermadec Islands",
        "South of Kermadec Islands",
        "North of Fiji Islands",
        "Fiji Islands region",
        "Fiji Islands",
        "Santa Cruz Islands region",
        "Santa Cruz Islands",
        "Vanuatu Islands region",
        "Vanuatu Islands",
        "New Caledonia",
        "Loyalty Islands",
        "Southeast of Loyalty Islands",
        "New Ireland region",
        "North of Solomon Islands",
        "New Britain region",
        "Bougainville - Solomon Islands region",
        "D'Entrecasteaux Islands region",
        "South of Solomon Islands",
        "Irian Jaya region",
        "Near north coast of Irian Jaya",
        "Ninigo Islands region",
        "Admiralty Islands region",
        "Near north coast of New Guinea",
        "Irian Jaya",
        "New Guinea",
        "Bismarck Sea",
        "Aru Islands region",
        "Near south coast of Irian Jaya",
        "Near south coast of New Guinea",
        "Eastern New Guinea region",
        "Arafura Sea",
        "Western Caroline Islands",
        "South of Mariana Islands",
        "Southeast of Honshu",
        "Bonin Islands region",
        "Volcano Islands region",
        "West of Mariana Islands",
        "Mariana Islands region",
        "Mariana Islands",
        "Kamchatka Peninsula",
        "Near east coast of Kamchatka Peninsula",
        "Off east coast of Kamchatka Peninsula",
        "Northwest of Kuril Islands",
        "Kuril Islands",
        "East of Kuril Islands",
        "Eastern Sea of Japan",
        "Hokkaido region",
        "Off southeast coast of Hokkaido",
        "Near west coast of eastern Honshu",
        "Eastern Honshu",
        "Near east coast of eastern Honshu",
        "Off east coast of Honshu",
        "Near south coast of eastern Honshu",
        "South Korea",
        "Western Honshu",
        "Near south coast of western Honshu",
        "Northwest of Ryukyu Islands",
        "Kyushu",
        "Shikoku",
        "Southeast of Shikoku",
        "Ryukyu Islands",
        "Southeast of Ryukyu Islands",
        "West of Bonin Islands",
        "Philippine Sea",
        "Near coast of southeastern China",
        "Taiwan region",
        "Taiwan",
        "Northeast of Taiwan",
        "Southwestern Ryukyu Islands",
        "Southeast of Taiwan",
        "Philippine Islands region",
        "Luzon",
        "Mindoro",
        "Samar",
        "Palawan",
        "Sulu Sea",
        "Panay",
        "Cebu",
        "Leyte",
        "Negros",
        "Sulu Archipelago",
        "Mindanao",
        "East of Philippine Islands",
        "Borneo",
        "Celebes Sea",
        "Talaud Islands",
        "North of Halmahera",
        "Minahassa Peninsula, Sulawesi",
        "Northern Molucca Sea",
        "Halmahera",
        "Sulawesi",
        "Southern Molucca Sea",
        "Ceram Sea",
        "Buru",
        "Seram",
        "Southwest of Sumatera",
        "Southern Sumatera",
        "Java Sea",
        "Sunda Strait",
        "Jawa",
        "Bali Sea",
        "Flores Sea",
        "Banda Sea",
        "Tanimbar Islands region",
        "South of Jawa",
        "Bali region",
        "South of Bali",
        "Sumbawa region",
        "Flores region",
        "Sumba region",
        "Savu Sea",
        "Timor region",
        "Timor Sea",
        "South of Sumbawa",
        "South of Sumba",
        "South of Timor",
        "Myanmar-India border region",
        "Myanmar-Bangladesh border region",
        "Myanmar",
        "Myanmar-China border region",
        "Near south coast of Myanmar",
        "Southeast Asia (REGION NOT IN USE)",
        "Hainan Island",
        "South China Sea",
        "Eastern Kashmir",
        "Kashmir-India border region",
        "Kashmir-Xizang border region",
        "Western Xizang-India border region",
        "Xizang",
        "Sichuan",
        "Northern India",
        "Nepal-India border region",
        "Nepal",
        "Sikkim",
        "Bhutan",
        "Eastern Xizang-India border region",
        "Southern India",
        "India-Bangladesh border region",
        "Bangladesh",
        "Northeastern India",
        "Yunnan",
        "Bay of Bengal",
        "Kyrgyzstan-Xinjiang border region",
        "Southern Xinjiang",
        "Gansu",
        "Western Nei Mongol",
        "Kashmir-Xinjiang border region",
        "Qinghai",
        "Southwestern Siberia",
        "Lake Baykal region",
        "East of Lake Baykal",
        "Eastern Kazakhstan",
        "Lake Issyk-Kul region",
        "Kazakhstan-Xinjiang border region",
        "Northern Xinjiang",
        "Tuva-Buryatia-Mongolia border region",
        "Mongolia",
        "Ural Mountains region",
        "Western Kazakhstan",
        "Eastern Caucasus",
        "Caspian Sea",
        "Northwestern Uzbekistan",
        "Turkmenistan",
        "Iran-Turkmenistan border region",
        "Turkmenistan-Afghanistan border region",
        "Turkey-Iran border region",
        "Iran-Armenia-Azerbaijan border region",
        "Northwestern Iran",
        "Iran-Iraq border region",
        "Western Iran",
        "Northern and central Iran",
        "Northwestern Afghanistan",
        "Southwestern Afghanistan",
        "Eastern Arabian Peninsula",
        "Persian Gulf",
        "Southern Iran",
        "Southwestern Pakistan",
        "Gulf of Oman",
        "Off coast of Pakistan",
        "Ukraine - Moldova - Southwestern Russia region",
        "Romania",
        "Bulgaria",
        "Black Sea",
        "Crimea region",
        "Western Caucasus",
        "Greece-Bulgaria border region",
        "Greece",
        "Aegean Sea",
        "Turkey",
        "Turkey-Georgia-Armenia border region",
        "Southern Greece",
        "Dodecanese Islands",
        "Crete",
        "Eastern Mediterranean Sea",
        "Cyprus region",
        "Dead Sea region",
        "Jordan - Syria region",
        "Iraq",
        "Portugal",
        "Spain",
        "Pyrenees",
        "Near south coast of France",
        "Corsica",
        "Central Italy",
        "Adriatic Sea",
        "Northwestern Balkan Peninsula",
        "West of Gibraltar",
        "Strait of Gibraltar",
        "Balearic Islands",
        "Western Mediterranean Sea",
        "Sardinia",
        "Tyrrhenian Sea",
        "Southern Italy",
        "Albania",
        "Greece-Albania border region",
        "Madeira Islands region",
        "Canary Islands region",
        "Morocco",
        "Northern Algeria",
        "Tunisia",
        "Sicily",
        "Ionian Sea",
        "Central Mediterranean Sea",
        "Near coast of Libya",
        "North Atlantic Ocean",
        "Northern Mid-Atlantic Ridge",
        "Azores Islands region",
        "Azores Islands",
        "Central Mid-Atlantic Ridge",
        "North of Ascension Island",
        "Ascension Island region",
        "South Atlantic Ocean",
        "Southern Mid-Atlantic Ridge",
        "Tristan da Cunha region",
        "Bouvet Island region",
        "Southwest of Africa",
        "Southeastern Atlantic Ocean",
        "Eastern Gulf of Aden",
        "Socotra region",
        "Arabian Sea",
        "Lakshadweep region",
        "Northeastern Somalia",
        "North Indian Ocean",
        "Carlsberg Ridge",
        "Maldive Islands region",
        "Laccadive Sea",
        "Sri Lanka",
        "South Indian Ocean",
        "Chagos Archipelago region",
        "Mauritius - Reunion region",
        "Southwest Indian Ridge",
        "Mid-Indian Ridge",
        "South of Africa",
        "Prince Edward Islands region",
        "Crozet Islands region",
        "Kerguelen Islands region",
        "Broken Ridge",
        "Southeast Indian Ridge",
        "Southern Kerguelen Plateau",
        "South of Australia",
        "Saskatchewan",
        "Manitoba",
        "Hudson Bay",
        "Ontario",
        "Hudson Strait region",
        "Northern Quebec",
        "Davis Strait",
        "Labrador",
        "Labrador Sea",
        "Southern Quebec",
        "Gaspe Peninsula",
        "Eastern Quebec",
        "Anticosti Island",
        "New Brunswick",
        "Nova Scotia",
        "Prince Edward Island",
        "Gulf of St. Lawrence",
        "Newfoundland",
        "Montana",
        "Eastern Idaho",
        "Hebgen Lake region, Montana",
        "Yellowstone region",
        "Wyoming",
        "North Dakota",
        "South Dakota",
        "Nebraska",
        "Minnesota",
        "Iowa",
        "Wisconsin",
        "Illinois",
        "Michigan",
        "Indiana",
        "Southern Ontario",
        "Ohio",
        "New York",
        "Pennsylvania",
        "Vermont - New Hampshire region",
        "Maine",
        "Southern New England",
        "Gulf of Maine",
        "Utah",
        "Colorado",
        "Kansas",
        "Iowa-Missouri border region",
        "Missouri-Kansas border region",
        "Missouri",
        "Missouri-Arkansas border region",
        "Missouri-Illinois border region",
        "New Madrid region, Missouri",
        "Cape Girardeau region, Missouri",
        "Southern Illinois",
        "Southern Indiana",
        "Kentucky",
        "West Virginia",
        "Virginia",
        "Chesapeake Bay region",
        "New Jersey",
        "Eastern Arizona",
        "New Mexico",
        "Northwestern Texas-Oklahoma border region",
        "Western Texas",
        "Oklahoma",
        "Central Texas",
        "Arkansas-Oklahoma border region",
        "Arkansas",
        "Louisiana-Texas border region",
        "Louisiana",
        "Mississippi",
        "Tennessee",
        "Alabama",
        "Western Florida",
        "Georgia",
        "Florida-Georgia border region",
        "South Carolina",
        "North Carolina",
        "Off east coast of United States",
        "Florida Peninsula",
        "Bahama Islands",
        "Eastern Arizona-Sonora border region",
        "New Mexico-Chihuahua border region",
        "Texas-Mexico border region",
        "Southern Texas",
        "Near coast of Texas",
        "Chihuahua",
        "Northern Mexico",
        "Central Mexico",
        "Jalisco",
        "Veracruz",
        "Gulf of Mexico",
        "Bay of Campeche",
        "Brazil",
        "Guyana",
        "Suriname",
        "French Guiana",
        "Eire",
        "United Kingdom",
        "North Sea",
        "Southern Norway",
        "Sweden",
        "Baltic Sea",
        "France",
        "Bay of Biscay",
        "The Netherlands",
        "Belgium",
        "Denmark",
        "Germany",
        "Switzerland",
        "Northern Italy",
        "Austria",
        "Czech and Slovak Republics",
        "Poland",
        "Hungary",
        "Northwest Africa (REGION NOT IN USE)",
        "Southern Algeria",
        "Libya",
        "Egypt",
        "Red Sea",
        "Western Arabian Peninsula",
        "Chad region",
        "Sudan",
        "Ethiopia",
        "Western Gulf of Aden",
        "Northwestern Somalia",
        "Off south coast of northwest Africa",
        "Cameroon",
        "Equatorial Guinea",
        "Central African Republic",
        "Gabon",
        "Congo",
        "Zaire",
        "Uganda",
        "Lake Victoria region",
        "Kenya",
        "Southern Somalia",
        "Lake Tanganyika region",
        "Tanzania",
        "Northwest of Madagascar",
        "Angola",
        "Zambia",
        "Malawi",
        "Namibia",
        "Botswana",
        "Zimbabwe",
        "Mozambique",
        "Mozambique Channel",
        "Madagascar",
        "South Africa",
        "Lesotho",
        "Swaziland",
        "Off coast of South Africa",
        "Northwest of Australia",
        "West of Australia",
        "Western Australia",
        "Northern Territory",
        "South Australia",
        "Gulf of Carpentaria",
        "Queensland",
        "Coral Sea",
        "Northwest of New Caledonia",
        "New Caledonia region",
        "Southwest of Australia",
        "Off south coast of Australia",
        "Near coast of South Australia",
        "New South Wales",
        "Victoria",
        "Near southeast coast of Australia",
        "Near east coast of Australia",
        "East of Australia",
        "Norfolk Island region",
        "Northwest of New Zealand",
        "Bass Strait",
        "Tasmania region",
        "Southeast of Australia",
        "North Pacific Ocean",
        "Hawaiian Islands region",
        "Hawaiian Islands",
        "Eastern Caroline Islands region",
        "Marshall Islands region",
        "Enewetak Atoll region",
        "Bikini Atoll region",
        "Gilbert Islands region",
        "Johnston Island region",
        "Line Islands region",
        "Palmyra Island region",
        "Kiritimati region",
        "Tuvalu region",
        "Phoenix Islands region",
        "Tokelau Islands region",
        "Northern Cook Islands",
        "Cook Islands region",
        "Society Islands region",
        "Tubuai Islands region",
        "Marquesas Islands region",
        "Tuamotu Archipelago region",
        "South Pacific Ocean",
        "Lomonosov Ridge",
        "Arctic Ocean",
        "Near north coast of Kalaallit Nunaat",
        "Eastern Kalaallit Nunaat",
        "Iceland region",
        "Iceland",
        "Jan Mayen Island region",
        "Greenland Sea",
        "North of Svalbard",
        "Norwegian Sea",
        "Svalbard region",
        "North of Franz Josef Land",
        "Franz Josef Land",
        "Northern Norway",
        "Barents Sea",
        "Novaya Zemlya",
        "Kara Sea",
        "Near coast of northwestern Siberia",
        "North of Severnaya Zemlya",
        "Severnaya Zemlya",
        "Near coast of northern Siberia",
        "East of Severnaya Zemlya",
        "Laptev Sea",
        "Southeastern Siberia",
        "Primorye-Northeastern China border region",
        "Northeastern China",
        "North Korea",
        "Sea of Japan",
        "Primorye",
        "Sakhalin Island",
        "Sea of Okhotsk",
        "Southeastern China",
        "Yellow Sea",
        "Off east coast of southeastern China",
        "North of New Siberian Islands",
        "New Siberian Islands",
        "Eastern Siberian Sea",
        "Near north coast of eastern Siberia",
        "Eastern Siberia",
        "Chukchi Sea",
        "Bering Strait",
        "St. Lawrence Island region",
        "Beaufort Sea",
        "Northern Alaska",
        "Northern Yukon Territory",
        "Queen Elizabeth Islands",
        "Northwest Territories",
        "Western Kalaallit Nunaat",
        "Baffin Bay",
        "Baffin Island region",
        "Southeast Central Pacific Ocean",
        "Southern East Pacific Rise",
        "Easter Island region",
        "West Chile Rise",
        "Juan Fernandez Islands region",
        "East of North Island",
        "Chatham Islands region",
        "South of Chatham Islands",
        "Pacific-Antarctic Ridge",
        "Southern Pacific Ocean",
        "East Central Pacific Ocean",
        "Central East Pacific Rise",
        "West of Galapagos Islands",
        "Galapagos Islands region",
        "Galapagos Islands",
        "Southwest of Galapagos Islands",
        "Southeast of Galapagos Islands",
        "South of Tasmania",
        "West of Macquarie Island",
        "Balleny Islands region",
        "Andaman Islands region",
        "Nicobar Islands region",
        "Off west coast of northern Sumatera",
        "Northern Sumatera",
        "Malay Peninsula",
        "Gulf of Thailand",
        "Southeastern Afghanistan",
        "Pakistan",
        "Southwestern Kashmir",
        "India-Pakistan border region",
        "Central Kazakhstan",
        "Southeastern Uzbekistan",
        "Tajikistan",
        "Kyrgyzstan",
        "Afghanistan-Tajikistan border region",
        "Hindu Kush region",
        "Tajikistan-Xinjiang border region",
        "Northwestern Kashmir",
        "Finland",
        "Norway-Murmansk border region",
        "Finland-Karelia border region",
        "Baltic States - Belarus - Northwestern Russia",
        "Northwestern Siberia",
        "Northern and central Siberia",
        "Victoria Land",
        "Ross Sea",
        "Antarctica",
        "Northern East Pacific Rise",
        "North of Honduras",
        "East of South Sandwich Islands",
        "Thailand",
        "Laos",
        "Cambodia",
        "Vietnam",
        "Gulf of Tongking",
        "Reykjanes Ridge",
        "Azores-Cape St. Vincent Ridge",
        "Owen Fracture Zone region",
        "Indian Ocean Triple Junction",
        "Western Indian-Antarctic Ridge",
        "Western Sahara",
        "Mauritania",
        "Mali",
        "Senegal - Gambia region",
        "Guinea region",
        "Sierra Leone",
        "Liberia region",
        "Cote d'Ivoire",
        "Burkina Faso",
        "Ghana",
        "Benin - Togo region",
        "Niger",
        "Nigeria",
        "Southeast of Easter Island",
        "Galapagos Triple Junction region",
    );
    return "UNKNOWN GEOGRAPHIC REGION" if ($grn < 1 || $grn > 757);
    return $regname[$grn - 1];
}

#
# $srname = sregion($srn)
#     returns seismic geographic name for a seismic region number
#  Input Arguments:
#     $srn - seismic region number
#  Return:
#     $srname - seismic region name
#
sub sregion {
    my $srn = shift;
    my @sregname = (
            "ALASKA - ALEUTIAN ARC", "E. ALASKA TO VANCOUVER ISLAND",
            "CALIFORNIA - NEVADA REGION", "BAJA CALIF. AND GULF OF CALIF.",
            "MEXICO - GUATEMALA AREA", "CENTRAL AMERICA", "CARIBBEAN LOOP",
            "ANDEAN SOUTH AMERICA", "EXTREME SOUTH AMERICA", "SOUTHERN ANTILLES",
            "NEW ZEALAND REGION", "KERMADEC - TONGA - SAMOA AREA",
            "FIJI ISLANDS AREA", "NEW HEBRIDES ISLANDS",
            "BISMARCK AND SOLOMON ISLANDS", "NEW GUINEA",
            "CAROLINE ISLANDS TO GUAM", "GUAM TO JAPAN",
            "JAPAN - KURILES - KAMCHATKA", "S.W. JAPAN AND RYUKYU ISLANDS",
            "TAIWAN", "PHILIPPINES", "BORNEO - CELEBES", "SUNDA ARC",
            "BURMA AND SOUTHEAST ASIA", "INDIA - TIBET - SZECHWAN - YUNAN",
            "SOUTHERN SINKIANG TO KANSU", "ALMA-ATA TO LAKE BAIKAL",
            "WESTERN ASIA", "MIDDLE EAST - CRIMEA - BALKANS",
            "WESTERN MEDITERRANEAN AREA", "ATLANTIC OCEAN", "INDIAN OCEAN",
            "EASTERN NORTH AMERICA", "EASTERN SOUTH AMERICA",
            "NORTHWESTERN EUROPE", "AFRICA", "AUSTRALIA", "PACIFIC BASIN",
            "ARCTIC ZONE", "EASTERN ASIA", "NE. ASIA,N. ALASKA TO GREENLAND",
            "S.E. AND ANTARCTIC PACIFIC", "GALAPAGOS AREA", "MACQUARIE LOOP",
            "ANDAMAN ISLANDS TO SUMATRA", "BALUCHISTAN", "HINDU KUSH AND PAMIR",
            "NORTHERN ASIA", "ANTARCTICA"
    );
    return "UNKNOWN SEISMIC REGION" if ($srn < 1 || $srn > 50);
    return $sregname[$srn - 1];
}

#
# $srn = gregtosreg($grn)
#     returns seismic region number for a geographic region number
#  Input Arguments:
#     $grn - geographic region number
#  Return:
#     $srn - seismic region number
#
sub gregtosreg {
    my $grn = shift;
    my ($ilo, $ihi);
    my @gstab = (0, 1, 18, 30, 47, 53, 72, 84, 102, 143,
                 147, 158, 169, 180, 183, 190, 196, 209, 211, 217,
                 231, 242, 248, 261, 273, 294, 302, 320, 326, 335,
                 357, 376, 402, 415, 438, 528, 532, 550, 588, 611,
                 633, 656, 667, 683, 693, 700, 703, 709, 713, 721,
                 727, 729);
    return -1 if ($grn < 1 || $grn > 757);
    if ($grn < 730) {
        ($ilo, $ihi) = bracket($grn, 51, \@gstab);
        return $ilo if ($grn < $gstab[$ihi]);
        return $ihi;
    }
    return  5 if ($grn == 730);
    return  7 if ($grn == 731);
    return 10 if ($grn == 732);
    return 25 if ($grn > 732 && $grn < 738);
    return 32 if ($grn > 737 && $grn < 740);
    return 33 if ($grn > 739 && $grn < 743);
    return 37 if ($grn > 742 && $grn < 756);
    return 43 if ($grn == 756);
    return 44 if ($grn == 757);
    return -1;
}

#
# $fe = read_FlinnEngdahl($fname)
#     reads Flinn-Engdahl regionalization file and stores it in a hash
#  Input Arguments:
#     $fname - pathname for Flinn-Engdahl regionalization file
#  Return:
#     $fe - hash pointer to Flinn-Engdahl regionalization
#  Flinn-Engdahl hash:
#     $fe->{error}  - set if error in reading file
#     $fe->{nlat}   - number of latitudes
#     @{$fe->{nl}}  - number of longitudes at each latitude(nlat)
#     @{$fe->{lon}} - array of longitudes(nlat x nl[i])
#     @{$fe->{grn}} - array of geographic region numbers(nlat x nl[i])
#
sub read_FlinnEngdahl {
    my $fname = shift;
    my ($fe, $nlat, $i, $j, $k, $m);
    my (@nl, @lon, @grn, @w);
    @nl = @lon = @grn = ();
    unless (open(IN, $fname)) {
        $fe->{error} = "No file: $fname";
        print "Cannot open $fname!\n";
        return $fe;
    }
    $k = $i = 0;
    while (<IN>) {
        next if (/^#/);
        if ($k == 0) {
#
#           number of latitudes
#
            ($nlat, undef) = split;
            $k++;
        }
        elsif ($k == 1) {
#
#           number of longitudes at each latitude
#
            @nl = split;
            $k++;
        }
        elsif ($k == 2) {
#
#           longitude ranges at each latitude
#
            ($m, @w) = split;
            for ($j = 0; $j < $nl[$i]; $j++) {
                $lon[$i][$j] = $w[$j];
            }
            $i++;
            if ($i == $nlat) {
                $i = 0;
                $k++;
            }
        }
        elsif ($k == 3) {
#
#           grn ranges at each latitude
#
            ($m, @w) = split;
            for ($j = 0; $j < $nl[$i]; $j++) {
                $grn[$i][$j] = $w[$j];
            }
            $i++;
            if ($i == $nlat) {
                $i = 0;
                $k++;
            }
        }
        else {
            next;
        }
    }
    close IN;
    $fe->{nlat} = $nlat;
    @{$fe->{nl}} = @nl;
    @{$fe->{lon}} = @lon;
    @{$fe->{grn}} = @grn;
    return $fe;
}

#
# $grn = gregnum($lat, $lon, $fe)
#     returns Flinn-Engdahl geographic region number for a lat, lon pair
#  Input Arguments:
#     $lat - latitude
#     $lon - longitude
#     $fe  - hash pointer to Flinn-Engdahl regionalization
#  Return:
#     $grn - Flinn-Engdahl geographic region number
#
sub gregnum {
    my ($lat, $lon, $fe) = @_;
    my ($ilat, $ilon, $ilo, $ihi, $n);
    $ilat = int(90. - $lat);
    $ilon = floor($lon);
    $n = $fe->{nl}[$ilat];
    return $fe->{grn}[$ilat][0] if ($n == 2);
    ($ilo, $ihi) = bracket($ilon, $n, $fe->{lon}[$ilat]);
    return $fe->{grn}[$ilat][$ilo] if ($ilon < $fe->{lon}[$ilat][$ihi]);
    return $fe->{grn}[$ilat][$ilo] if ($ilon == 180);
    return $fe->{grn}[$ilat][$ihi];
}

1;
