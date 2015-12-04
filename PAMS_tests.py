import struct
import time
import PAMS_func
import os
import re
import numpy as np

def Isodat_File_Parser(fileName):
    '''Reads in a .did file (Isodat acquisition file), returns the raw voltages for
    ref gas, analysis gas, and the isodat-calculated d13C and d18O of the Acquisition'''

    f=open(fileName,'rb')
    try:
        buff = f.read()
    finally:
        f.close()

    #1. Getting raw voltage data
    #Searching for the start of the raw voltages
    start=buff.find('CIntensityData')
    keys=[]
    voltRef_raw=[]
    voltSam_raw=[]
    #Slightly elegant method: pulling out based on spacing after 'CIntensityData' ascii key
    # TODO: make this more flexible
    # TODO: Catch errors if wrong voltage sequence found
    startPreVolt=start+4164 #observed location of 'pre' cycle on ref gas side
    voltRef_raw.append(struct.unpack('10d',buff[startPreVolt:(startPreVolt+10*8)]))
    for i in range(10):
        startRefVolt=start+52+i*200 #observed location of ref gas voltage cycles
        voltRef_raw.append(struct.unpack('10d',buff[startRefVolt:(startRefVolt+10*8)]))

        startSamVolt=start+2196+i*196 #observed location of sample gas voltage cycles
        voltSam_raw.append(struct.unpack('10d',buff[startSamVolt:(startSamVolt+10*8)]))


    # 3.2 whether or not method is a CO2_multiply or a *_start
    firstAcq = False
    startMethod = buff.find('CDualInletBlockData')
    methodBlock = buff[startMethod-120:startMethod-20].decode('utf-16')
    # if 'CO2_multiply_16V' in methodBlock:
    #     firstAcq = False
    # if 'Contin_Start' in methodBlock:
    #     lastAcq = True
    #     #TODO: make this a more robust test
    if 'methane' not in methodBlock:
        print('WARNING: Possibly incorrect acquisition method')

    # 3.3 sample name
    # Find start of block with sample name
    startName = buff.find('CSeqLineIndexData')
    # Rough guess of range where analysis name is, accounting for a large variation in length
    nameBlock = buff[startName+200:startName+400].decode('utf-16')
    #Exact name based on locations of unicode strings directly before and after
    analysisName = nameBlock[(nameBlock.find('Pressadjust')+19):(nameBlock.find('Identifier')-2)]
    # Encode as ascii for consistency
    analysisName = analysisName.encode('ascii', 'replace')

    # 3.4 background values
    #find start of block with background values
    # startBackground = buff.find('CISLScriptMessageData')
    # stopBackground = buff.find('CMeasurmentErrors')
    # #Note incorrect spelling of 'measurement' is intentional
    # backgroundBlock = buff[startBackground+80:stopBackground].decode('utf-16')

    # 3.5 Date and Time
    # Find the start of Ctime block
    startTime = buff.find('CTimeObject')+49
    # Pull out the time_t time based on startTime location, seconds since the epoch (Jan 1st, 1970), GMT
    time_t_time = struct.unpack('i',buff[startTime:startTime+4])[0]
    # time_str = time.strftime('%m/%d/%Y %H:%M:%S', time.localtime(time_t_time))
    date_str = time.strftime('%m/%d/%Y', time.localtime(time_t_time))
    time_str = time.strftime('%H:%M', time.localtime(time_t_time))



    return voltRef_raw, voltSam_raw, analysisName, date_str, time_str

analyses = []

while True:
    analyses.append(PAMS_func.MCI())
    while True:
        acqName = raw_input('Drag an acq file for sample ')
        if len(acqName) == 0:
            break

        # acqNum=re.findall('[0-9]{4}',acqName.split('/')[-1])[0] #finds the acquision number from the file name, no matter how long the path nor what it contains
        acqNum=re.findall('[0-9]{4}',os.path.basename(acqName))[0] #finds the acquision number from the file name, no matter how long the path nor what it contains
        acqNum=int(acqNum)
        acqName=acqName.strip()


        voltRef_raw, voltSam_raw, analysisName, date_str, time_str, deuteriumMeasurement, time_c = PAMS_func.Isodat_File_Parser(acqName)
        print('Name = {0}, date is {1}, time is {2}').format(analysisName,date_str,time_str)
        print('Deuterium measurement is {0}?').format(str(deuteriumMeasurement))
        print(voltRef_raw[-1])
        print(voltSam_raw[-1])
        if deuteriumMeasurement:
            analyses[-1].acqs_D.append(PAMS_func.ACQUISITION_D(acqNum))
            analyses[-1].acqs_D[-1].voltSam_raw = voltSam_raw
            analyses[-1].acqs_D[-1].voltRef_raw = voltRef_raw
            analyses[-1].acqs_D[-1].date = date_str
            analyses[-1].acqs_D[-1].time = time_str
            analyses[-1].acqs_D[-1].time_c = time_c

        else:
            analyses[-1].acqs_full.append(PAMS_func.ACQUISITION_FULL(acqNum))
            analyses[-1].acqs_full[-1].voltSam_raw = voltSam_raw
            analyses[-1].acqs_full[-1].voltRef_raw = voltRef_raw
            analyses[-1].acqs_full[-1].date = date_str
            analyses[-1].acqs_full[-1].time = time_str
            analyses[-1].acqs_full[-1].time_c = time_c

    newSample = raw_input('do another sample? (y/n)').lower()
    if len(newSample) == 0:
        break
