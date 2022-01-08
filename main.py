from conc_parser import process_data
import time

begin = time.time()
process_data(['data/LADDERROOM-A.plot','data/LADDERROOM-B.plot','data/LOOPSEAL-ROOM.plot'],'data/Temp.plot','data/volume.dat')
end = time.time()
print("Time taken is ", end - begin)