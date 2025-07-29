import cdsapi

dataset = "reanalysis-era5-pressure-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "geopotential",
        "specific_humidity",
        "temperature",
        "u_component_of_wind",
        "v_component_of_wind",
        "vertical_velocity"
    ],
    "year": ["2005"],
    "month": ["08"],
    "day": [
        "08", "09", "10",
        "11", "12", "13",
        "14",
    ],
    "time": [
        "00:00", "06:00", "12:00",
        "18:00"
    ],
    "pressure_level": [
        "10", "20", "30",
        "50", "70", "100",
        "125", "150", "175",
        "200", "225", "250",
        "300", "350", "400",
        "450", "500", "550",
        "600", "650", "700",
        "750", "775", "800",
        "825", "850", "875",
        "900", "925", "950",
        "975", "1000"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [-17.5, -60, -42.5, -30]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download("system-20050808_ERA5.nc")