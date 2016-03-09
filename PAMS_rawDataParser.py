'''Program to process raw isodat intensities for clumped CO2 measurements, and turn them into meaningful isotope ratios'''


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import PAMS_func
import os


analyses = []
print('Welcome to the methane clumped isotope importer/exporter')
while True:
    print('Please select a task:')
    taskChoice = raw_input(' (I)mport from raw isodat files \n Import from (P)AMS file \n (E)xport to PAMS and FLATLIST \n (C)alculate stretching corr and hg line \n (Q)uit \n ').upper()

    if taskChoice == 'I':
        print('This script turns raw isodat files into a FLATLIST ')
        # list containing instances of CI class, each with their own acq objects
        # list containing numbers of acquisitions already imported
        imported=[]
        modeChoice=raw_input('(a)utomatic mode or (m)anual mode? ').lower()
        if modeChoice == 'm':
            print('Manual mode selected')
            while True:
                sampleName = raw_input('Name of new sample, or press RETURN to stop: ')
                if len(sampleName) == 0:
                    break
                analyses.append(PAMS_func.MCI())
                analyses[-1].name=sampleName
                analyses[-1].num=acqNum
                while True:
                    acqName = raw_input('Drag an acq file for sample ' + analyses[-1].name +', or press RETURN to stop: ')
                    acqName=acqName.strip()

                    if len(acqName) == 0:
                        break
                    acqName = acqName.strip('"')
                    acqName = os.path.abspath(acqName)
                    # acqNum=re.findall('[0-9]{4}',acqName.split('/')[-1])[0] #finds the acquision number from the file name, no matter how long the path nor what it contains
                    acqNum=re.findall('[0-9]{4}',os.path.basename(acqName))[0] #finds the acquision number from the file name, no matter how long the path nor what it contains
                    acqNum=int(acqNum)
                    if acqNum in imported:
                        print('You already imported this file')
                    else:
                        voltRef_raw,voltSam_raw,rawSampleName, date_str, time_str, deuteriumMeasurement, time_c = PAMS_func.Isodat_File_Parser(acqName)
                        if sampleName != rawSampleName :
                            print('Sample name: ' + analyses[-1].name + ' does not match name in file: ' + rawSampleName + ' ')
                            nameErrorChoice = raw_input('Are you sure you want to include this acquisition (y/n)? ')
                            if nameErrorChoice.lower() == 'n':
                                continue

                        imported.append(acqNum)
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

                        print('Acquisition '+str(acqNum)+ ' successfully imported.')
                # Catching situation where sample name used in error, so no acqs imported
                if len(analyses[-1].acqs) == 0:
                    print('No acqs imported for sample ' + analyses[-1].name)
                    print('Deleting sample ' + analyses[-1].name +'...')
                    del analyses[-1]

        elif modeChoice == 'a':
            print('Automatic mode selected')
            while True:
                acqFolder = raw_input('Drag a folder containing all acquisitions: ').strip()
                acqFolder = acqFolder.strip('"')
                acqFolder = os.path.abspath(acqFolder)
                if os.path.exists(acqFolder):
                    break
                else:
                    print('Invalid folder, please try again ')

            acqList = [i for i in os.listdir(acqFolder) if i.endswith('.did')]
            if len(acqList) == 0:
                print("No acquisiton files ('.did') found in folder ")
                quit()
            acqList=sorted(acqList)
            startNum = raw_input('Number of first acquisition to be processed, or (f)irst: ').lower()
            stopNum = raw_input('Number of last acquisition to be processed, or (l)ast: ').lower()
            #convert to int first to remove any leading zeros
            if startNum in ['f', '0', '']:
                startNum = int(acqList[0].rstrip('.did').replace('Acquisition-',''))
            else:
                startNum = int(startNum)
            if stopNum in ['l', '-1', 'end', '']:
                stopNum = int(acqList[-1].rstrip('.did').replace('Acquisition-',''))
            else:
                stopNum = int(stopNum)

            #now, adding the proper number of leading zeros in order to get an exact match
            startName = 'Acquisition-' + (4-len(str(startNum)))*'0' + str(startNum) + '.did'
            stopName = 'Acquisition-' + (4-len(str(stopNum)))*'0' + str(stopNum) + '.did'

            startNumIndex = acqList.index(startName)
            stopNumIndex = acqList.index(stopName)

            firstAcq = False
            for i in range(startNumIndex,stopNumIndex+1):
                acqName = acqFolder +'/'+ acqList[i]
                # Catches files with a size less than 130 kb and skips them
                if os.path.getsize(acqName) < 130000:
                    print('Skipping acq num ' + str(acqList[i]) + ' because file too small')
                    continue
                # Finds the acquision number from the file name, no matter how long the path nor what it contains
                acqNum=re.findall('[0-9]{4}',acqName.split('/')[-1])[0]
                acqNum=int(acqNum)
                # Actually processing the file
                print('Importing acq num ' + str(acqNum) + ' ')
                voltRef_raw, voltSam_raw, rawSampleName, date_str, time_str, deuteriumMeasurement, time_c = PAMS_func.Isodat_File_Parser(acqName)
                # Creates a new sample if acquisition no acquisition already

                # Catches acqs where enough gas did not make it to the bellows in and skips them
                # if voltSam_raw[-1][0] < 15000:
                #     print('Skipping acq ' + str(acqList[i]) + ' from sample ' + rawSampleName + ' because voltage too low on mass 44: ' + str(voltSam_raw[-1][0]))
                #     continue
                # if first acquision of folder, new sample
                if len(analyses) == 0:
                    analyses.append(PAMS_func.MCI())
                    analyses[-1].name = rawSampleName
                    analyses[-1].num = acqNum
                    analyses[-1].date = date_str
                    lastSampleName = rawSampleName
                    firstAcq = False
                    print('Found new sample, with name: ' + rawSampleName)
                # if had been measuring clumped species, and now measuring D, new samples
                if deuteriumMeasurement and len(analyses[-1].acqs_D)>0:
                    tempAcq = PAMS_func.ACQUISITION_D(acqNum)
                    tempAcq.voltSam_raw = voltSam_raw
                    tempAcq.voltRef_raw = voltRef_raw
                    d17_temp = tempAcq.d17_D
                    d17_last = analyses[-1].acqs_D[-1].d17_D
                    if abs(d17_temp-d17_last) > 8 and analyses[-1].name.lower() != rawSampleName.lower():
                        firstAcq = True
                elif not deuteriumMeasurement and len(analyses[-1].acqs_full)>0:
                    tempAcq = PAMS_func.ACQUISITION_FULL(acqNum)
                    tempAcq.voltSam_raw = voltSam_raw
                    tempAcq.voltRef_raw = voltRef_raw
                    d17_temp = tempAcq.d17_full
                    d17_last = analyses[-1].acqs_full[-1].d17_full
                    if abs(d17_temp-d17_last) > 1 and analyses[-1].name.lower() != rawSampleName.lower():
                        firstAcq = True

                if firstAcq:
                    analyses.append(PAMS_func.MCI())
                    analyses[-1].name = rawSampleName
                    analyses[-1].num = acqNum
                    analyses[-1].date = date_str
                    firstAcq = False
                    print('Found new sample, with name: ' + rawSampleName)
                # Catching case where name in file does not match current acq name

                if not firstAcq and rawSampleName.lower() not in [lastSampleName.lower(), analyses[-1].name.lower()] :
                    print('Old sample name: ' + analyses[-1].name + ' does not match name in this file: ' + rawSampleName + ' ')
                    print('Last d17 measurement: {0:.3f}, this d17 measurement: {1:.3f}'.format(d17_last, d17_temp))
                    nameErrorChoice = raw_input('(s)kip acquisition, (i)nclude it, or make a (n)ew sample from it? ')
                    if nameErrorChoice.lower() == 's':
                        print('Skipping acquisition ')
                        continue
                    elif nameErrorChoice.lower() == 'n':
                        print('Making a new sample with name: ' + rawSampleName)
                        analyses.append(PAMS_func.MCI())
                        analyses[-1].name = rawSampleName
                        analyses[-1].num = acqNum
                        analyses[-1].date = date_str
                        firstAcq = False
                    else:
                        print('Including acquisition ')
                # Catching case where just d17s don't match up
                # elif not firstAcq:
                #     if abs(d17_temp-d17_last) > (8*deuteriumMeasurement + 1*(not deuteriumMeasurement)):
                #         print('d17 mismatch detected in sample {0} '.format(rawSampleName))
                #         print('Last d17 measurement: {0:.3f}, this d17 measurement: {1:.3f}'.format(d17_last, d17_temp))
                #         nameErrorChoice = raw_input('(s)kip acquisition, (i)nclude it, or make a (n)ew sample from it? ')
                #         if nameErrorChoice.lower() == 's':
                #             print('Skipping acquisition ')
                #             continue
                #         elif nameErrorChoice.lower() == 'n':
                #             # print('Making a new sample with name: ' + rawSampleName)
                #             # analyses.append(PAMS_func.MCI())
                #             # analyses[-1].name = rawSampleName
                #             # analyses[-1].num = acqNum
                #             # analyses[-1].date = date_str
                #             firstAcq = True
                #         else:
                #             print('Including acquisition ')
                #
                # if firstAcq:
                #     analyses.append(PAMS_func.MCI())
                #     analyses[-1].name = rawSampleName
                #     analyses[-1].num = acqNum
                #     analyses[-1].date = date_str
                #     firstAcq = False
                #     print('Found new sample, with name: ' + rawSampleName)


                # if no errors caught above, actually add acquisition to analyses
                if deuteriumMeasurement:
                    analyses[-1].acqs_D.append(PAMS_func.ACQUISITION_D(acqNum))
                    analyses[-1].acqs_D[-1].voltSam_raw = voltSam_raw
                    analyses[-1].acqs_D[-1].voltRef_raw = voltRef_raw
                    analyses[-1].acqs_D[-1].date = date_str
                    analyses[-1].acqs_D[-1].time = time_str
                    analyses[-1].acqs_D[-1].name = rawSampleName
                    analyses[-1].acqs_D[-1].time_c = time_c

                else:
                    analyses[-1].acqs_full.append(PAMS_func.ACQUISITION_FULL(acqNum))
                    analyses[-1].acqs_full[-1].voltSam_raw = voltSam_raw
                    analyses[-1].acqs_full[-1].voltRef_raw = voltRef_raw
                    analyses[-1].acqs_full[-1].date = date_str
                    analyses[-1].acqs_full[-1].time = time_str
                    analyses[-1].acqs_full[-1].name = rawSampleName
                    analyses[-1].acqs_full[-1].time_c = time_c

                lastSampleName = rawSampleName

                print('Acquisition '+str(acqNum)+ ' successfully imported.')

        else :
            print('Not a valid mode choice')
            print('Goodbye')

        if len(analyses) != 0:
            print('Acquisition imports complete')
            print(str(len(analyses)) + ' analyses were imported')
            print('Cleaning up analyses...')
            analyses=PAMS_func.PAMS_cleaner(analyses)
            print('Checking analyses types...')
            while not PAMS_func.Sample_type_checker(analyses):
                print('Some analyses types need to be assigned ')
                typeGetterMode = raw_input('Get analyses types (a)utomatically or (m)anually? ').lower()
                if typeGetterMode == 'm':
                    analyses = PAMS_func.Get_types_manual(analyses)
                else:
                    analyses = PAMS_func.Get_types_auto(analyses)

            # exportNowChoice = raw_input('Export analyses now (y/n)? ').lower()
            # if exportNowChoice == 'y':
            #     PAMS_func.ExportSequence(analyses)

    #
    if taskChoice == 'P':
        print('Importing analyses from one or more auto-formatted PAMS files')
        while True:
            filePath = raw_input('Drag in a file or q to stop importing: ').strip().lower()
            if (filePath == 'q' or len(filePath) == 0):
                break
            if '.csv' not in filePath:
                print('file must be a .csv format')
                continue
            filePath = filePath.strip('"')
            filePath = os.path.abspath(filePath)
            print('Importing file... ')
            newAnalyses = PAMS_func.PAMS_importer(filePath)
            print('{0} new analyses imported from file'.format(len(newAnalyses)))
            analyses += newAnalyses
            print('Cleaning up analyses...')
            analyses=PAMS_func.PAMS_cleaner(analyses)
            print('Checking analyses types...')
            while not PAMS_func.Sample_type_checker(analyses):
                print('Some analyses types need to be assigned ')
                typeGetterMode = raw_input('Get analyses types (a)utomatically or (m)anually? ').lower()
                if typeGetterMode == 'm':
                    analyses = PAMS_func.Get_types_manual(analyses)
                else:
                    analyses = PAMS_func.Get_types_auto(analyses)

            print('Successfully added to dataset')
        print('{0} total analyses imported'.format(len(analyses)))

    if taskChoice == 'E':
        if len(analyses) == 0:
            print('Nothing to export!')
            break
        else:
            while not PAMS_func.Sample_type_checker(analyses):
                print('Some analyses types need to be assigned ')
                typeGetterMode = raw_input('Get analyses types (a)utomatically or (m)anually? ').lower()
                if typeGetterMode == 'm':
                    analyses = PAMS_func.Get_types_manual(analyses)
                else:
                    analyses = PAMS_func.Get_types_auto(analyses)

            PAMS_func.ExportSequence(analyses)
        #
    if taskChoice == 'C':
        print('Determining correct stretching values')
        PAMS_func.MCI_stretching_determination(analyses, showFigures = True)
        print('Processing data in heated gas ref frame')
        PAMS_func.MCI_hg_data_corrector(analyses)
    if taskChoice == 'Q':
        print('Goodbye! ')
        break
    else:
        pass


























#
#
# # For drag-and-drop capacity on windows systems
# # Ask for file name otherwise
# if (len(sys.argv) < 2):
#     fileName = raw_input('Name of the file to parse ? ')
# else:
#     fileName = sys.argv[1]
#
# FILE = open(FILENAME, 'r')
# FILENAME = FILENAME.rstrip('.txt')
#
#
# samples=CIDS_func.CIDS_parser(filePath)
#
# samples =CIDS_func.CIDS_cleaner(samples)
# samples=CIDS_func.D47_calculations(samples)
#
# pblChoice = raw_input('Would you like to do a pressure baseline correction? (y/n) ')
#
# if pblChoice == 'y':
#     print 'Pressure baseline correction selected'
#     print ''
