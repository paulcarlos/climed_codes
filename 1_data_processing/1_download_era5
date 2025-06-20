import cdsapi

pnm = "" # add folder to save
dataset = "reanalysis-era5-single-levels"
#dataset = "reanalysis-era5-land"
varnm = "2m_dewpoint_temperature" # "total_precipitation" "2m_temperature" "2m_dewpoint_temperature"
varfo = "dp" # "tp" "t2m" "dp"
yrv = ["2000","2001","2002","2003","2004","2005"] # add years
mov = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

client = cdsapi.Client() # needs separate setup

for j in yrv:
    for i in mov:
        request = {
            "product_type": ["reanalysis"], # silent for ERA5-Land
            "variable": [varnm],
            "year": j,
            "month": i,
            "day": [
                "01", "02", "03",
                "04", "05", "06",
                "07", "08", "09",
                "10", "11", "12",
                "13", "14", "15",
                "16", "17", "18",
                "19", "20", "21",
                "22", "23", "24",
                "25", "26", "27",
                "28", "29", "30",
                "31"
            ],
            "time": [
                "00:00", "01:00", "02:00",
                "03:00", "04:00", "05:00",
                "06:00", "07:00", "08:00",
                "09:00", "10:00", "11:00",
                "12:00", "13:00", "14:00",
                "15:00", "16:00", "17:00",
                "18:00", "19:00", "20:00",
                "21:00", "22:00", "23:00"
            ],
            "data_format": "grib",
            "download_format": "unarchived"
        }
        fnm = pnm + varfo + "/" + j + "/" + varfo + "_hourly_era5_global_" + j + "-" + i + ".grib"
        client.retrieve(dataset, request, fnm)
