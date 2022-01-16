from conc_parser import process_data

process_data(
    ["data/LADDERROOM-A.plot", "data/LADDERROOM-B.plot", "data/LOOPSEAL-ROOM.plot"],
    "data/Temp.plot",
    "data/volume.dat",
)
