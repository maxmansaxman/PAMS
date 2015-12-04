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
    taskChoice = raw_input(' (I)mport from raw isodat files \n Import from (P)AMS file \n (E)xport to PAMS and FLATLIST \n (A)nalyze data \n (Q)uit \n ').upper()

    if taskChoice == 'I':
        print('This script turns raw isodat files into a FLATLIST ')
        analyses = []
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
                        voltRef_raw,voltSam_raw,rawSampleName, date_str, time_str, deuteriumMeasurement = PAMS_func.Isodat_File_Parser(acqName)
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
                        else:
                            analyses[-1].acqs_full.append(PAMS_func.ACQUISITION_FULL(acqNum))
                            analyses[-1].acqs_full[-1].voltSam_raw = voltSam_raw
                            analyses[-1].acqs_full[-1].voltRef_raw = voltRef_raw
                            analyses[-1].acqs_full[-1].date = date_str
                            analyses[-1].acqs_full[-1].time = time_str
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
            startNum = raw_input('Number of first acquisition to be processed: ')
            stopNum = raw_input('Number of last acquisition to be processed: ')
            #convert to int first to remove any leading zeros
            startNum = int(startNum)
            stopNum = int(stopNum)

            #now, adding the proper number of leading zeros in order to get an exact match
            startName = 'Acquisition-' + (4-len(str(startNum)))*'0' + str(startNum) + '.did'
            stopName = 'Acquisition-' + (4-len(str(stopNum)))*'0' + str(stopNum) + '.did'

            startNumIndex = acqList.index(startName)
            stopNumIndex = acqList.index(stopName)

            firstAcq = False
            measuringClumped = False
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
                    firstAcq = True
                # if had been measuring clumped species, and now measuring D, new samples
                if measuringClumped and deuteriumMeasurement:
                    firstAcq = True

                if firstAcq:
                    analyses.append(PAMS_func.MCI())
                    analyses[-1].name = rawSampleName
                    analyses[-1].num = acqNum
                    analyses[-1].date = date_str
                    firstAcq = False
                    print('Found new sample, with name: ' + rawSampleName)

                # Catching case where name in file does not match current acq name
                if analyses[-1].name.lower() != rawSampleName.lower() :
                    print('Sample name: ' + analyses[-1].name + ' does not match name in file: ' + rawSampleName + ' ')
                    nameErrorChoice = raw_input('(s)kip acquisition, (i)nclude it, or make a (n)ew sample from it? ')
                    if nameErrorChoice.lower() == 's':
                        print('Skipping acquisition ')
                        continue
                    elif nameErrorChoice.lower() == 'n':
                        print('Making a new sample with name: ' + rawSampleName)
                        analyses.append(PAMS_func.MCI())
                        analyses[-1].name = rawSampleName
                        analyses[-1].num = acqNum
                        analyses[-1].date = date
                        firstAcq = False
                    else:
                        print('Including acquisition ')


                # if no errors caught above, actually add acquisition to analyses
                if deuteriumMeasurement:
                    analyses[-1].acqs_D.append(PAMS_func.ACQUISITION_D(acqNum))
                    analyses[-1].acqs_D[-1].voltSam_raw = voltSam_raw
                    analyses[-1].acqs_D[-1].voltRef_raw = voltRef_raw
                    analyses[-1].acqs_D[-1].date = date_str
                    analyses[-1].acqs_D[-1].time = time_str
                    analyses[-1].acqs_D[-1].name = rawSampleName
                    analyses[-1].acqs_D[-1].time_c = time_c
                    measuringClumped = False

                else:
                    analyses[-1].acqs_full.append(PAMS_func.ACQUISITION_FULL(acqNum))
                    analyses[-1].acqs_full[-1].voltSam_raw = voltSam_raw
                    analyses[-1].acqs_full[-1].voltRef_raw = voltRef_raw
                    analyses[-1].acqs_full[-1].date = date_str
                    analyses[-1].acqs_full[-1].time = time_str
                    analyses[-1].acqs_full[-1].name = rawSampleName
                    analyses[-1].acqs_full[-1].time_c = time_c
                    measuringClumped = True

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

            exportNowChoice = raw_input('Export analyses now (y/n)? ').lower()
            if exportNowChoice == 'y':
                PAMS_func.ExportSequence(analyses)

    #
    # if taskChoice == 'C':
    #     print('Importing analyses from one or more auto-formatted CIDS files')
    #     while True:
    #         filePath = raw_input('Drag in a file or q to stop importing: ').strip().lower()
    #         if (filePath == 'q' or len(filePath) == 0):
    #             break
    #         if '.csv' not in filePath:
    #             print('file must be a .csv format')
    #             continue
    #         filePath = filePath.strip('"')
    #         filePath = os.path.abspath(filePath)
    #         print('Importing file... ')
    #         newAnalyses = CIDS_func.CIDS_importer(filePath)
    #         print('{0} new analyses imported from file'.format(len(newAnalyses)))
    #         analyses += newAnalyses
    #         print('Checking analyses types...')
    #         while not CIDS_func.Sample_type_checker(analyses):
    #             print('Some analyses types need to be assigned ')
    #             typeGetterMode = raw_input('Get analyses types (a)utomatically or (m)anually? ').lower()
    #             if typeGetterMode == 'm':
    #                 analyses = CIDS_func.Get_types_manual(analyses)
    #             else:
    #                 analyses = CIDS_func.Get_types_auto(analyses)
    #
    #         print('Successfully added to dataset')
    #     print('{0} total analyses imported'.format(len(analyses)))

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
    # if taskChoice == 'P':
    #     print('Processing data in all relevant reference frames')
    #     print('Processing data in caltech ref frame')
    #     CIDS_func.CI_CRF_data_corrector(analyses)
    #     print('Processing data in absolute ref frame, with the Daeron method')
    #     CIDS_func.Daeron_data_processer(analyses)
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