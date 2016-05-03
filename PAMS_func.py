''' Stores all the functions and classes used by various clumped isotope scripts'''


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import struct
import time
import xlcor47_modified
import string
import datetime
import pandas as pd
from scipy.optimize import root
from scipy.optimize import brentq
from scipy.interpolate import interp1d

class MCI_VALUE(object):
    '''subclass defining how important isotopic ratios are calculated'''
    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if len(instance.voltRef)>3:
            if self.name in ['d17_full', 'd18']:
                return np.around(MCI_calculation_full(instance, self.name),3)
            elif self.name in ['d17_D']:
                return np.around(MCI_calculation_D(instance, self.name),3)

    def __set__(self, obj, value):
        raise AttributeError('Cannot change MCI calculation scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete MCI value')

class MCI_AVERAGE(object):
    '''subclass defining how isotopic ratios are averaged'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if len(instance.acqs_full)>=1:
            if self.name in ['d17_full','d17_full_sterr','d17_D','d17_D_sterr','d18','d18_sterr']:
                return np.around(MCI_averages(instance, self.name),3)
        else:
            raise ValueError('Sample has no acquisitions to average')

    def __set__(self, obj, value):
        raise AttributeError('Cannot change MCI averaging scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete MCI average')

class MCI_CALCULATED_VALUE(object):
    '''subclass defining how isotopic ratios are averaged'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if self.name in ['d13C', 'd13C_sterr', 'dD', 'dD_sterr', 'D18_raw', 'D18_sterr']:
            return np.around(MCI_bulk_comp(instance, self.name),3)

        else:
            raise ValueError('Sample D47_raw is out of range')

    def __set__(self, obj, value):
        raise AttributeError('Cannot change CI correction scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete CI corretion scheme')
#
class MCI_APPLIED_VALUE(object):
    '''subclass defining how applied vales are set in all acqs'''

    def __init__(self, name, initval):
        self.name = name
        self.val = initval

    def __get__(self,instance,cls):
        if self.name == 'stretch_17':
            return instance.stretch_17_holder
        elif self.name == 'stretch_18':
            return instance.stretch_18_holder
        else:
            raise ValueError('Not a valid value')

    def __set__(self, instance, value):
        if self.name in ['stretch_17', 'stretch_18']:
            self.val = value
            MCI_apply_to_acqs(instance, self.name, self.val)

    def __delete__(self, instance):
        raise AttributeError('Cannot delete individual stretching value')

class MCI_CORRECTED_VALUE(object):
    '''subclass defining how isotopic ratios are averaged'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if self.name in ['D18_hg']:
            return np.around(MCI_hg_corrector(instance, self.name),3)
        elif self.name in ['D18_stoch']:
            return np.around(instance.D18_hg + 2.981,3)
        else:
            raise ValueError('Sample D47_raw is out of range')

    def __set__(self, obj, value):
        raise AttributeError('Cannot change CI correction scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete CI corretion scheme')


class MCI_VOLTAGE(object):
    '''subclass defining how raw voltages are corrected'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if self.name in ['voltSam', 'voltRef']:
            return np.around(MCI_background_correction(instance, self.name),3)

        else:
            raise ValueError('Not a valid value')

    def __set__(self, obj, value):
        raise AttributeError('Cannot change individual voltage correction scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete individual voltage correction scheme')

class MCI(object):
    "A class for all the attributes of a single clumped isotope measurement"

    def __init__(self):
        self.name=''
        self.acqs_full=[]
        self.acqs_D=[]
        self.date = ''
        self.type=''
        self.num=np.nan
        self.user=''

        self.d13C_ref = -42.88
        self.dD_ref = -175.5

        self.D18 = np.nan
        self.D48_error = np.nan

        self.adduct_full = [0, 0]
        self.adduct_D = [0, 0]
        self.adduct_18 = [0, 0]

        self.stretch_18_holder = 3.700725
        self.stretch_17_holder =  0.053245
        self.frag = 0.789

        self.hg_slope = 0
        self.hg_intercept = 0

        self.background_full = [0,0,0]
        self.background_D = [0,0,0]

    stretch_17 = MCI_APPLIED_VALUE('stretch_17', 0.053245)
    stretch_18 = MCI_APPLIED_VALUE('stretch_18', 3.700725)

    d17_full = MCI_AVERAGE('d17_full')
    d18 = MCI_AVERAGE('d18')
    d17_D = MCI_AVERAGE('d17_D')

    D18_raw = MCI_CALCULATED_VALUE('D18_raw')
    d13C = MCI_CALCULATED_VALUE('d13C')
    dD = MCI_CALCULATED_VALUE('dD')

    d17_full_sterr = MCI_AVERAGE('d17_full_sterr')
    d18_sterr = MCI_AVERAGE('d18_sterr')
    d17_D_sterr = MCI_AVERAGE('d17_D_sterr')

    D18_sterr = MCI_CALCULATED_VALUE('D18_sterr')
    d13C_sterr = MCI_CALCULATED_VALUE('d13C_sterr')
    dD_sterr = MCI_CALCULATED_VALUE('dD_sterr')

    D18_hg = MCI_CORRECTED_VALUE('D18_hg')
    D18_stoch = MCI_CORRECTED_VALUE('D18_stoch') # Value in stochastic ref frame
    # hg_slope = MCI_UNIVERSAL_VALUE('hg_slope')
    # hg_intercept = MCI_UNIVERSAL_VALUE('hg_intercept')

class ACQUISITION_FULL(object):
    "A class for all the attributes of a single clumped isotope acquision"

    def __init__(self,acqNum):
        self.acqNum=acqNum
        self.voltSam_raw=[]
        self.voltRef_raw=[]
        self.background=[0, 0, 0]
        self.date=''
        self.time=''
        self.adduct_17 = [ 0, 0]
        self.adduct_18 = [ 0, 0]
        self.stretch_18 = 3.700725
        self.type = 'full'
        self.name = ''
        self.time_c = 0

    d17_full=MCI_VALUE('d17_full')
    d18=MCI_VALUE('d18')
    voltSam = MCI_VOLTAGE('voltSam')
    voltRef = MCI_VOLTAGE('voltRef')

class ACQUISITION_D(object):
    "A class for all the attributes of a single clumped isotope acquision"
    def __init__(self,acqNum):
        self.acqNum=acqNum
        self.voltSam_raw=[]
        self.voltRef_raw=[]
        self.background=[0, 0, 0]
        self.date=''
        self.time=''
        self.adduct_17 = [0, 0]
        self.stretch_17= 0.053245
        self.type = 'D'
        self.time_c = 0
        self.name = ''

    d17_D=MCI_VALUE('d17_D')
    voltSam = MCI_VOLTAGE('voltSam')
    voltRef = MCI_VOLTAGE('voltRef')

def Is_It_MethaneShutdown(fileName):
    ''' Reads in a .did file and decides whether the method called was methaneshutdown'''
    f = open(fileName, 'rb')
    try:
        buff = f.read()
    finally: f.close()
    startMethod = buff.find('CDualInletBlockData')
    methodBlock = buff[startMethod-120:startMethod-20].decode('utf-16')

    if 'methaneshutdown' in methodBlock:
        return(True)
    else:
        return(False)


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

    # 3.2 whether or not method is a methane method
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
    # function to filter out non-ascii chars that are tacked on to end
    # and, for that matter, anything else following these
    anf = ''
    for i in analysisName:
        if i not in string.printable:
            break
        anf += i
    # Encode as ascii for consistency
    analysisName = anf.encode('ascii', 'ignore')

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

    # 3.6 Deciding if d13C or dD measurement
    # converting to arrays
    voltRef_raw = np.asarray(voltRef_raw)
    voltSam_raw = np.asarray(voltSam_raw)

    if np.mean(voltRef_raw[:,7]) < 5e3:
        deuteriumMeasurement = True
    else:
        deuteriumMeasurement = False

    return voltRef_raw, voltSam_raw, analysisName, date_str, time_str, deuteriumMeasurement, time_t_time

def PAMS_cleaner(analyses):
    '''function for cleaning up a parsed CIDS file and alerting to any corrupted analyses'''
    # Cleaning up the parsed file
    # because of the way the parser is written, an extra empty analysis is added to the end of the list
    if not analyses[-1].acqs_full:
        del analyses[-1]
    # checking for analyses with not enough acqs, alerting if there are some
    temp_full=8
    lowAcqs= [k for k in analyses if len(k.acqs_full)<temp_full]
    if not lowAcqs:
        print 'All analyses have at least %d clumped acquistions' % temp_full

    else:
        print 'analyses with too few clumped acqs are:'
        print '\n'.join(['Sample '+ str(k.num)+ ': '+ k.name +' has ' + str(len(k.acqs_full)) + ' acquisitions' for k in lowAcqs])
        for i in lowAcqs:
            if len(i.acqs_full) <=1:
                analyses.remove(i)
    temp_D=3
    lowAcqs_D= [k for k in analyses if len(k.acqs_D)<temp_D]
    if not lowAcqs_D:
        print 'All analyses have at least %d D acquistions' % temp_D

    else:
        print 'analyses with too few D acqs are:'
        print '\n'.join(['Sample '+ str(k.num)+ ': '+ k.name +' has ' + str(len(k.acqs_D)) + ' acquisitions' for k in lowAcqs_D])
        for i in lowAcqs_D:
            if len(i.acqs_D) <=1:
                analyses.remove(i)
    # Setting background corrections
    checkBackgroundsChoice = raw_input('Set background corrections now? (y/n) ').lower()
    if checkBackgroundsChoice == 'y':
        print('Checking background lines')
        for i in analyses:
            backVals = [i.background_D, i.background_full]
            backTypes = ['dD', 'full']
            for j in range(2):
                if backVals[j] == [0,0,0]:
                    print('Sample {0} from {1} needs {2} background vals '.format(i.name, i.date, backTypes[j]))
                    a = raw_input('mass 16: ')
                    a = float(a)
                    b = raw_input('mass 17: ')
                    b = float(b)
                    c = 0.0
                    if j ==1:
                        c = raw_input('mass 18: ')
                        c = float(c)

                    # temporarily store values
                    backVals[j] = [a, b, c]
            # now, setting values
            [i.background_D, i.background_full] = backVals

    # Setting adduct lines
    checkAdductLinesChoice = raw_input('Set adduct lines now? (y/n) ').lower()
    if checkAdductLinesChoice == 'y':
        print('Checking adduct lines')
        for i in analyses:
            adductVals = [i.adduct_D, i.adduct_full, i.adduct_18]
            adductTypes = ['D', 'full 17', 'clumped 18']
            for j in range(3):
                if adductVals[j] == [0,0]:
                    print('Sample {0} from {1} needs a {2} line '.format(i.name, i.date, adductTypes[j]))
                    b = raw_input('1st-order coefficient: ')
                    b = float(b)
                    a = raw_input('2nd-order coefficient: ')
                    a = float(a)
                    # temporarily store values
                    adductVals[j] = [a, b]
            # now, setting values
            [i.adduct_D, i.adduct_full, i.adduct_18] = adductVals
    # Regardless of status, set all acq-level adduct lines to sample adduct lines
    # AND, set all acq-level backgrounds to sample level backgrounds
    # NOTE: This assumes that backgrounds are uniform for all acqs of the same type in each sample.
    for i in analyses:
        for k in i.acqs_D:
            k.adduct_17 = i.adduct_D
            k.background = i.background_D
        for l in i.acqs_full:
            [l.adduct_17, l.adduct_18] = [i.adduct_full, i.adduct_18]
            l.background = i.background_full
    # Setting stretching corrections
    changeStretchingCorrections = raw_input('Change stretching corrections now? (y/n) ' ).lower()
    if changeStretchingCorrections == 'y':
        for i in analyses:
            print('Now considering sample {0} from {1}:  (press q to quit)'.format(i.name, i.date))
            print('dD stretching correction is currently {0}'.format(i.stretch_17))
            changeIt = raw_input('Input new one or press enter to keep: ').lower()
            if changeIt == 'q':
                break
            elif len(changeIt) > 0:
                i.stretch_17 = float(changeIt)

            print('d18 stretching correction is currently {0}'.format(i.stretch_18))
            changeIt2 = raw_input('Input new one or press enter to keep: ').lower()
            if len(changeIt2) > 0:
                i.stretch_18 = float(changeIt2)

    print 'All analyses are cleaned, and voltages converted to arrays'

    return analyses

def FlatList_exporter(analyses,fileName, displayProgress = False):
    '''Exports a CSV file with the essential calculated MCI info'''

    export=open(fileName + '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    wrt.writerow(['User','date','Type','Sample ID','spec #\'s', 'full acqs', 'dD acqs', 'd13C (vpdb)','d13C_sterr','d2H (vsmow)',
    'd2H_sterr','d18','d18_sterr','D18 (v. wg)','D18_sterr', 'hg_slope', 'hg_intercept', 'D18_stochastic'])
    counter = 0
    progressBar = np.linspace(0, len(analyses), 20)
    for item in analyses:
        wrt.writerow([item.user, item.date, item.type, item.name, item.num, len(item.acqs_full), len(item.acqs_D), item.d13C, item.d13C_sterr, item.dD,
        item.dD_sterr,item.d18,item.d18_sterr,item.D18_raw,item.D18_sterr, item.hg_slope, item.hg_intercept, item.D18_stoch ])
        counter += 1
        if displayProgress:
            if float(counter) in np.floor(progressBar):
                print(str(int(float(counter)/len(analyses)*100)) + '% done')

    export.close()
    return

def PAMS_exporter(analyses, fileName, displayProgress = False):
    '''Exports a CSV file that is the same format as a traditional CIDS file.
    Python CIDS files are also importable, to load all important sample info
    back into the program'''

    export=open(fileName + '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    counter = 0
    progressBar = np.linspace(0, len(analyses), 20)
    for analysis in analyses:
        wrt.writerow(['__NewSample__'])
        for acq in analysis.acqs_D:
            wrt.writerow(['__NewAcq-D__'])
            wrt.writerow(['__AcqVoltage__','',]+np.shape(acq.voltSam_raw)[1]*['',]+acq.voltRef_raw[0,:].tolist())
            for i in range(len(acq.voltSam_raw)):
                wrt.writerow(['',''] + acq.voltSam_raw[i,:].tolist() + acq.voltRef_raw[i+1,:].tolist())
            wrt.writerow([])
            wrt.writerow(['','acq num','date','name','type', 'mass16_bckgrnd', 'mass17_bckgrnd','adductLine (2nd-ordr)', 'adducLine (1st-ordr)','17stretchCorrection','d17_D (wg)'])
            wrt.writerow(['__AcqInfo__',acq.acqNum, acq.date, analysis.name, acq.type, acq.background[0], acq.background[1],acq.adduct_17[0], acq.adduct_17[1], acq.stretch_17, acq.d17_D])

        wrt.writerow([])

        for acq in analysis.acqs_full:
            wrt.writerow(['__NewAcq-Full__'])
            wrt.writerow(['__AcqVoltage__','',]+np.shape(acq.voltSam_raw)[1]*['',]+acq.voltRef_raw[0,:].tolist())
            for i in range(len(acq.voltSam_raw)):
                wrt.writerow(['',''] + acq.voltSam_raw[i,:].tolist() + acq.voltRef_raw[i+1,:].tolist())
            wrt.writerow([])
            wrt.writerow(['','acq num','date','name','type', 'mass16_bckgrnd', 'mass17_bckgrnd', 'mass18_bckgrnd', '17adductLine (2nd-ordr)',
            '17adductLine (1st-ordr)', '18adductLine (2nd-ordr)', '18adductLine (1st-ordr)','18stretchCorrection','d17_full (wg)', 'd18 (wg)'])
            wrt.writerow(['__AcqInfo__',acq.acqNum, acq.date, analysis.name, acq.type, acq.background[0], acq.background[1], acq.background[2], acq.adduct_17[0], acq.adduct_17[1],
            acq.adduct_18[0], acq.adduct_18[1],acq.stretch_18, acq.d17_full, acq.d18])

        wrt.writerow([])
        wrt.writerow(['D analyses', 'd17_D','','','full analyses','d17_full','d18'])
        for i in range(max([len(analysis.acqs_D),len(analysis.acqs_full)])):
            try:
                d17_D = analysis.acqs_D[i].d17_D
            except(IndexError):
                d17_D = ''
            try:
                d17_full = analysis.acqs_full[i].d17_full
                d18 = analysis.acqs_full[i].d18
            except(IndexError):
                d17_full = ''
                d18 = ''
            wrt.writerow(['',d17_D, '', '', '' , d17_full, d18])

        wrt.writerow([])
        wrt.writerow(['','User','date','Type','Sample ID','acq #\'s','adduct_D (2nd-ordr)',
        'adduct_D (1st-ordr)','adduct_full (2nd-ordr)','adduct_full (1st-ordr)','adduct_18 (2nd-ordr)',
        'adduct_18 (1st-ordr)','D_stretchCorrection','18_stretchCorrection','Fragmenation','dD_bckgrnd_mass16', 'dD_bckgrnd_mass17', 'full_bckgrnd_mass16',
        'full_bckgrnd_mass17', 'full_bckgrnd_mass18','d13C (vpdb)','d13C_sterr','dD (vsmow)','dD_sterr','D18_raw','D18_sterr', 'hg_slope', 'hg_intercept', 'D18_hg'])
        wrt.writerow(['__SampleSummary__',analysis.user, analysis.date, analysis.type, analysis.name, analysis.num, analysis.adduct_D[0], analysis.adduct_D[1],
        analysis.adduct_full[0], analysis.adduct_full[1], analysis.adduct_18[0], analysis.adduct_18[1], analysis.stretch_17, analysis.stretch_18, analysis.frag,
        analysis.background_D[0], analysis.background_D[1], analysis.background_full[0], analysis.background_full[1], analysis.background_full[2], analysis.d13C,
        analysis.d13C_sterr, analysis.dD, analysis.dD_sterr, analysis.D18_raw, analysis.D18_sterr, analysis.hg_slope, analysis.hg_intercept, analysis.D18_hg])
        wrt.writerow(31*['---',])

        counter +=1
        if displayProgress:
            if float(counter) in np.floor(progressBar):
                print(str(int(float(counter)/len(analyses)*100)) + '% done')

    export.close()
    return

def PAMS_importer(filePath, displayProgress = False):
    '''Imports voltage and analysis data that were exported using the PAMS_exporter function'''
    fileImport = open(filePath,'rU')
    fileReader = csv.reader(fileImport, dialect='excel')
    analyses = []

    cycles_in_acqs = 10

    for line in fileReader:
        if '__NewSample__' in line:
            analyses.append(MCI())
            continue
        if '__NewAcq-D__' in line:
            analyses[-1].acqs_D.append(ACQUISITION_D(0))
            importingDeuterium = True
            continue
        if '__NewAcq-Full__' in line:
            analyses[-1].acqs_full.append(ACQUISITION_FULL(0))
            importingDeuterium = False
            continue
        if '__AcqVoltage__' in line:
            voltIndex = line.index('__AcqVoltage__')
            cycle = 0
            voltRefs = [float(vr) for vr in line[voltIndex+12:voltIndex+22]]
            if importingDeuterium:
                analyses[-1].acqs_D[-1].voltRef_raw.append(voltRefs)
            else:
                analyses[-1].acqs_full[-1].voltRef_raw.append(voltRefs)
            continue
        if cycle < (cycles_in_acqs):
            voltSams=[float(vs) for vs in line[voltIndex+2:voltIndex+12]]
            voltRefs = [float(vr) for vr in line[voltIndex+12:voltIndex+22]]
            if importingDeuterium:
                analyses[-1].acqs_D[-1].voltSam_raw.append(voltSams)
                analyses[-1].acqs_D[-1].voltRef_raw.append(voltRefs)
            else:
                analyses[-1].acqs_full[-1].voltSam_raw.append(voltSams)
                analyses[-1].acqs_full[-1].voltRef_raw.append(voltRefs)
            cycle += 1
            continue
        if '__AcqInfo__' in line:
            acqIndex = line.index('__AcqInfo__')
            measType = line[acqIndex+4]
            if measType == 'full':
                analyses[-1].acqs_full[-1].acqNum = float(line[acqIndex + 1])
                analyses[-1].acqs_full[-1].date = line[acqIndex + 2]
                analyses[-1].acqs_full[-1].voltRef_raw = np.asarray(analyses[-1].acqs_full[-1].voltRef_raw)
                analyses[-1].acqs_full[-1].voltSam_raw = np.asarray(analyses[-1].acqs_full[-1].voltSam_raw)
                # analyses[-1].acqs_full[-1].background = [float(i) for i in line[acqIndex+5: acqIndex + 8]]

            elif measType == 'D':
                analyses[-1].acqs_D[-1].acqNum = float(line[acqIndex + 1])
                analyses[-1].acqs_D[-1].date = line[acqIndex + 2]
                analyses[-1].acqs_D[-1].voltRef_raw = np.asarray(analyses[-1].acqs_D[-1].voltRef_raw)
                analyses[-1].acqs_D[-1].voltSam_raw = np.asarray(analyses[-1].acqs_D[-1].voltSam_raw)
                # analyses[-1].acqs_D[-1].background = [float(i) for i in line[acqIndex+5: acqIndex + 7]].append(0.0)

            continue
        if '__SampleSummary__' in line:
            summaryIndex = line.index('__SampleSummary__')
            [analyses[-1].user, analyses[-1].date, analyses[-1].type, analyses[-1].name] = line[summaryIndex + 1:summaryIndex+5]
            analyses[-1].num = int(line[summaryIndex + 5])
            [analyses[-1].adduct_D[0],analyses[-1].adduct_D[1],analyses[-1].adduct_full[0],analyses[-1].adduct_full[1], analyses[-1].adduct_18[0],analyses[-1].adduct_18[1]] = [float(i) for i in line[summaryIndex + 6:summaryIndex+12]]
            [analyses[-1].stretch_17,analyses[-1].stretch_18, analyses[-1].frag] = [float(i) for i in line[summaryIndex+12:summaryIndex+15]]
            [analyses[-1].background_D[0], analyses[-1].background_D[1]] = [float(i) for i in line[summaryIndex+15:summaryIndex+17]]
            analyses[-1].background_full = [float(i) for i in line[summaryIndex+17:summaryIndex+20]]
            for i in analyses[-1].acqs_D:
                i.adduct_17 = analyses[-1].adduct_D
                i.stretch_17 = analyses[-1].stretch_17
            for j in analyses[-1].acqs_full:
                j.adduct_17 = analyses[-1].adduct_full
                j.adduct_18 = analyses[-1].adduct_18
                j.stretch_18 = analyses[-1].stretch_18

    fileImport.close()
    return analyses

def MCI_stretching_determination(analyses, showFigures = True):
    ''' Fuction to determine the stretching corrections with the +1 std and possibly apply it '''
    # accepted values of dD and D18  of +1 std
    std1_true_old = [226.92, 230.88] # Before April 2016
    std1_true_new = [320.531, 327.17] # After April 14th, 2016

    # list of all std+1 standards
    std1s = [i for i in analyses if i.type =='std1']
    dates = [i.date for i in std1s]

    dates_pd = pd.to_datetime(dates)
    # Between April 1st and April 14th, 2016, two differed standards were run that were both called +1 std
    # On April 14th, the switchover to the new one was finished
    if '2016-04' in dates_pd:
        print('WARNING: Value of +1 std changed in April 2016. \nThis dataset contains +1 stds from this month. \nAssuming these are newer values. \nBeware!')

    s17s = []
    s17s_upper= []
    s17s_lower = []
    s18s = []
    s18s_upper = []
    s18s_lower = []
    # Determining stretching values for each +1 standard
    for j in std1s:
        # Check for old or new std +1
        if dates_pd[i] < np.datetime64('2016-04-01'):
            std1_true = std1_true_old
        else:
            std1_true = std1_true_new
        # extra arg is this specific to this std+1
        extraArgs = (j,std1_true,)
        j.stretch_17 = brentq(MCI_std1_stretch_17_fn, 0.0001, 0.1, args = extraArgs, xtol = 1e-12, disp = False)
        j.stretch_18 = brentq(MCI_std1_stretch_18_fn, 0.0001, 10, args = extraArgs, xtol = 1e-12, disp = False)
        s17s.append(j.stretch_17)
        s18s.append(j.stretch_18)
        s17s_upper.append(brentq(MCI_std1_stretch_17_upper_fn, 0.0001, 0.1, args = extraArgs, xtol = 1e-12, disp = False))
        s17s_lower.append(brentq(MCI_std1_stretch_17_lower_fn, 0.0001, 0.1, args = extraArgs, xtol = 1e-12, disp = False))
        s18s_upper.append(brentq(MCI_std1_stretch_18_upper_fn, 0.0001, 10, args = extraArgs, xtol = 1e-12, disp = False))
        s18s_lower.append(brentq(MCI_std1_stretch_18_lower_fn, 0.0001, 10, args = extraArgs, xtol = 1e-12, disp = False))

    # Calculating mean errors
    s17sterr = []
    s18sterr = []
    for i in range(len(s17s_upper)):
        s17sterr.append(np.mean([np.abs(s17s_upper[i]-s17s[i]), np.abs(s17s_lower[i] - s17s[i])]))
        s18sterr.append(np.mean([np.abs(s18s[i]-s18s_upper[i]), np.abs(s18s[i] - s18s_lower[i])]))


    # Converting s17s into pandas data series
    s17s = pd.Series(np.asarray(s17s), index = dates_pd)
    s18s = pd.Series(np.asarray(s18s), index = dates_pd)
    s17sterr = pd.Series(np.asarray(s17sterr), index = dates_pd)
    s18sterr = pd.Series(np.asarray(s18sterr), index = dates_pd)
    # reindexing to add in all na values from all analyses
    dates_all = [i.date for i in analyses]
    # pandas date index
    dates_pd_all = pd.to_datetime(dates_all)
    # ctime series
    ctimes_all = []
    ctimes = []
    for i in dates_all:
        try:
            ctimes_all.append(time.mktime(time.strptime(i, '%m/%d/%Y')))
        except(ValueError):
            ctimes_all.append(time.mktime(time.strptime(i, '%m/%d/%y')))
        if i in dates:
            try:
                ctimes.append(time.mktime(time.strptime(i, '%m/%d/%Y')))
            except(ValueError):
                ctimes.append(time.mktime(time.strptime(i, '%m/%d/%y')))

    # ctimes = [time.mktime(time.strptime(i, '%m/%d/%Y')) for i in dates_all]
    s17s = s17s.reindex(dates_pd_all)
    s18s = s18s.reindex(dates_pd_all)
    s17sterr = s17sterr.reindex(dates_pd_all)
    s18sterr = s18sterr.reindex(dates_pd_all)

    if showFigures:
        plt.close('all')
        plt.ion()
        fig1, ax1 = plt.subplots()
        s17s.plot(ax = ax1, style = 'bo')
        s17s.plot(ax = ax1, style = 'bo-', yerr = s17sterr)
        ax1.set_xlabel('date')
        ax1.set_ylabel('stretch_17 of Std+1')
        ax1.set_xlim(s17s.index.min(), s17s.index.max())

        fig0, ax0 = plt.subplots()
        s18s.plot(ax = ax0, style = 'rd')
        s18s.plot(ax = ax0, style = 'rd-', yerr = s18sterr)
        ax0.set_xlabel('date')
        ax0.set_ylabel('stretch_18 of Std+1')
        ax0.set_xlim(s18s.index.min(), s18s.index.max())
        #
        # fig2, ax2 = plt.subplots()
        # ax2.plot(s17s, s18s, 'bo')
        # ax2.set_xlabel('stretch_17 of Std+1')
        # ax2.set_ylabel('stretch_18 of Std+1')

    stretchChoice = raw_input('Apply (s)ingle stretching correction, linear (r)egression, linear (i)nterpolation, or (f)orward fill? ').lower()
    if stretchChoice == 's':
        print('Mean stretch_17 is {0:.5f} '.format(s17s.mean()))
        stretchChoice2 = raw_input('Use (m)ean, or type another: ')
        if stretchChoice2 in ['m','', 'M']:
            chosen17 = s17s.mean()
        else:
            chosen17 = float(stretchChoice2)
        print('Mean stretch_18 is {0:.5f} '.format(s18s.mean()))
        stretchChoice3 = raw_input('Use (m)ean, or type another: ')
        if stretchChoice3 in ['m','', 'M']:
            chosen18 = s18s.mean()
        else:
            chosen18 = float(stretchChoice3)
        # Now, doing fill. Not changing std1s
        for i in analyses:
            if i.type != 'std1':
                i.stretch_17 = chosen17
                i.stretch_18 = chosen18

    elif stretchChoice == 'r':
        stretchFn_17 = np.polyfit(np.array(ctimes), s17s.dropna().values, 1)
        stretchFn_18 = np.polyfit(np.array(ctimes), s18s.dropna().values, 1)
        for i in range(len(analyses)):
            # if analyses[i].type != 'std1':
                # i.stretch_17 = time.mktime(time.strptime(i.date, '%m/%d/%y'))*stretchFn_17[0] + stretchFn_17[1]
                # i.stretch_18 = time.mktime(time.strptime(i.date, '%m/%d/%y'))*stretchFn_18[0] + stretchFn_18[1]
            analyses[i].stretch_17 = ctimes_all[i]*stretchFn_17[0] + stretchFn_17[1]
            analyses[i].stretch_18 = ctimes_all[i]*stretchFn_18[0] + stretchFn_18[1]

    elif stretchChoice == 'i':
        s17s.interpolate(method = 'time', inplace = True)
        s18s.interpolate(method = 'time', inplace = True)
        # now, back fill for first values
        s17s.fillna(method = 'bfill', inplace = True)
        s18s.fillna(method = 'bfill', inplace = True)
        # stretchFn_17 = interp1d(ctimes, s17s, kind = 'linear')
        # stretchFn_18 = interp1d(ctimes, s18s, kind = 'linear')
        # Now, fill into data set
        for i in range(len(analyses)):
            analyses[i].stretch_17 = s17s.ix[i]
            analyses[i].stretch_18 = s18s.ix[i]
            # analyses[i].stretch_17 = stretchFn_17(time.mktime(time.strptime(i.date, '%m/%d/%y')))
            # analyses[i].stretch_18 = stretchFn_18(time.mktime(time.strptime(i.date, '%m/%d/%y')))
    # elif stretchChoice == 'm':
    #     windowSize = raw_input('Size of moving window (2-10)? ').strip()
    #     windowSize = int(windowSize)
    #     s17s = pd.rolling_mean(s17s, windowSize, min_periods = 1)
    #     s18s = pd.rolling_mean(s18s, windowSize, min_periods = 1)
    #     # Now, fill into data set
    #     for i in range(len(analyses)):
    #         analyses[i].stretch_17 = s17s.ix[i]
    #         analyses[i].stretch_18 = s18s.ix[i]

    elif stretchChoice =='f':
        s17s.fillna(method = 'ffill', inplace = True)
        s18s.fillna(method = 'ffill', inplace = True)
        # Now, fill into data set
        for i in range(len(analyses)):
            analyses[i].stretch_17 = s17s.ix[i]
            analyses[i].stretch_18 = s18s.ix[i]

    # if showFigures:
    #     # Now, plot std results:
    #     hgs = [i for i in analyses if (i.type == 'hg')]
    #     d18_hgs = np.asarray([i.d18 for i in hgs])
    #     D18_raw_hgs = np.asarray([i.D18_raw for i in hgs])
    #     d18_sterr_hgs = np.asarray([i.d18_sterr for i in hgs])
    #     D18_sterr_hgs = np.asarray([i.D18_sterr for i in hgs])
    #
    #     stds_D13 = [i for i in analyses if (i.type == 'stdD13')]
    #     stds_1= [i for i in analyses if (i.type == 'std1')]
    #     fig1, ax1 = plt.subplots()
    #     ax1.errorbar([i.d18 for i in stds_D13],[i.D18_hg for i in stds_D13], xerr = [i.d18_sterr for i in stds_D13], yerr = [i.D18_sterr for i in stds_D13], fmt = 'o', label = '+D+13C')
    #     ax1.errorbar([i.d18 for i in stds_1],[i.D18_hg for i in stds_1], xerr = [i.d18_sterr for i in stds_1], yerr = [i.D18_sterr for i in stds_1], fmt = 's', label = '+1')
    #     ax1.errorbar([i.d18 for i in hgs],[i.D18_hg for i in hgs], xerr = [i.d18_sterr for i in hgs], yerr = [i.D18_sterr for i in hgs], fmt = '^', label = 'hg')
    #     ax1.set_xlabel(ur'$\mathrm{}\delta^{18} \/ (\u2030)}$')
    #     ax1.set_ylabel(ur'$\mathrm{\Delta_{18, hg_corr} \/ (\u2030)}$')
    #     ax1.legend()


    return


def MCI_std1_stretch_17_fn(stretchGuess, *extraArgs):
    this_std1, std1_true = extraArgs
    this_std1.stretch_17 = stretchGuess
    return(this_std1.d17_D - std1_true[0])

def MCI_std1_stretch_17_upper_fn(stretchGuess, *extraArgs):
    this_std1, std1_true = extraArgs
    this_std1.stretch_17 = stretchGuess
    return((this_std1.d17_D + this_std1.d17_D_sterr) - std1_true[0])

def MCI_std1_stretch_17_lower_fn(stretchGuess, *extraArgs):
    this_std1, std1_true = extraArgs
    this_std1.stretch_17 = stretchGuess
    return((this_std1.d17_D - this_std1.d17_D_sterr) - std1_true[0])

def MCI_std1_stretch_18_fn(stretchGuess, *extraArgs):
    this_std1, std1_true = extraArgs
    this_std1.stretch_18 = stretchGuess
    return(this_std1.d18 - std1_true[1])

def MCI_std1_stretch_18_upper_fn(stretchGuess, *extraArgs):
    this_std1, std1_true = extraArgs
    this_std1.stretch_18 = stretchGuess
    return((this_std1.d18 + this_std1.d18_sterr) - std1_true[1])

def MCI_std1_stretch_18_lower_fn(stretchGuess, *extraArgs):
    this_std1, std1_true = extraArgs
    this_std1.stretch_18 = stretchGuess
    return((this_std1.d18 - this_std1.d18_sterr) - std1_true[1])




def MCI_hg_data_corrector(analyses, showFigures = True):
    ''' Function to determine heated gas line, plot it, and apply it to all samples '''

    # First, make a list of all hgs

    hgs = [i for i in analyses if (i.type == 'hg')]
    d18_hgs = np.asarray([i.d18 for i in hgs])
    D18_raw_hgs = np.asarray([i.D18_raw for i in hgs])
    d18_sterr_hgs = np.asarray([i.d18_sterr for i in hgs])
    D18_sterr_hgs = np.asarray([i.D18_sterr for i in hgs])

    # Correcting for offset between HGs and known wg value
    D18_raw_hgs += 1.97
    # Now, make hg D18 line, using a York regression
    # global hg_slope, hg_intercept
    if len(hgs) <= 2:
        print('Not enough heated gases for a meaningful hg line! ')
        return
    hg_slope_temp, hg_intercept_temp, r_18, sm, sb, xc, yc, ct = lsqcubic(d18_hgs, D18_raw_hgs,d18_sterr_hgs, D18_sterr_hgs)

    D18_raw_hgs_model = d18_hgs*hg_slope_temp+hg_intercept_temp
    if showFigures:
        plt.ion()
        plt.figure(0)
        plt.errorbar(d18_hgs, D18_raw_hgs, xerr = d18_sterr_hgs, yerr = D18_sterr_hgs, fmt = 'bo')
        plt.plot(d18_hgs, D18_raw_hgs_model, '-')
        plt.ylabel(ur'$\Delta_{18} \/ (\u2030)$')
        plt.xlabel(ur'$\delta^{18} \/ (\u2030)$')


        # Calculation of the D18_hg values are done automatically given its definition, based on the global hg slope and int
        # Plot the stds and heated gases
        stds_D13 = [i for i in analyses if (i.type == 'stdD13')]
        stds_1= [i for i in analyses if (i.type == 'std1')]
        fig1, ax1 = plt.subplots()
        ax1.errorbar([i.d18 for i in stds_D13],[i.D18_raw for i in stds_D13], xerr = [i.d18_sterr for i in stds_D13], yerr = [i.D18_sterr for i in stds_D13], fmt = 'o', label = '+D+13C')
        ax1.errorbar([i.d18 for i in stds_1],[i.D18_raw for i in stds_1], xerr = [i.d18_sterr for i in stds_1], yerr = [i.D18_sterr for i in stds_1], fmt = 's', label = '+1')
        ax1.errorbar([i.d18 for i in hgs],[i.D18_raw for i in hgs], xerr = [i.d18_sterr for i in hgs], yerr = [i.D18_sterr for i in hgs], fmt = '^', label = 'hg')
        ax1.set_xlabel(ur'$\mathrm{}\delta^{18} \/ (\u2030)}$')
        ax1.set_ylabel(ur'$\mathrm{\Delta_{18, raw} \/ (\u2030)}$')
        ax1.legend()

    time.sleep(1)
    hgChoice = raw_input(' Apply this heated gas correction to all analyses? \n Slope: {0:.3f}, Intercept: {1:.3f} (y/n) '.format(hg_slope_temp, hg_intercept_temp)).lower()
    if hgChoice == 'y':
        print('Adding to all analyses... ')
        for i in analyses:
            i.hg_slope = hg_slope_temp
            i.hg_intercept = hg_intercept_temp
    return


def MCI_hg_slope_finder(hgs):
    # hgs = [i for i in analyses if (i.type == 'hg' and not i.D48_excess)]
    d18_hgs = np.asarray([i.d18 for i in hgs])
    D18_raw_hgs = np.asarray([i.D18_raw for i in hgs])
    d18_sterr_hgs = np.asarray([i.d18_sterr for i in hgs])
    D18_sterr_hgs = np.asarray([i.D18_sterr for i in hgs])
    # Now, make hg D18 line, using a York regression
    global hg_slope, hg_intercept
    hg_slope, hg_intercept, r_18, sm, sb, xc, yc, ct = lsqcubic(d18_hgs, D18_raw_hgs,d18_sterr_hgs, D18_sterr_hgs)

    return hg_slope

def MCI_apply_to_acqs(analysis, objName, objValue):
    ''' Function to apply a shared trait to all relevant acqs'''
    if objName  == 'stretch_17':
        analysis.stretch_17_holder = objValue
        for i in analysis.acqs_D:
            i.stretch_17 = objValue
    elif objName == 'stretch_18':
        analysis.stretch_18_holder = objValue
        for i in analysis.acqs_full:
            i.stretch_18 = objValue
    return


def MCI_hg_corrector(analysis, objName):
    '''Function to apply the heated gas correction in the Caltech Ref Frame'''

    try:
        D18_hg_corrected = analysis.D18_raw - (analysis.d18*analysis.hg_slope + analysis.hg_intercept)
        # D47_stretching = D47_hg_corrected *(-0.8453)/hg_intercept
        # D47_CRF = D47_stretching + acid_digestion_correction
    except NameError:
        return(np.NaN)

    return(D18_hg_corrected)

# def MCI_hg_values(instance, objName):
#     '''Placeholder for hg slope and int when used in CI classes'''
#     try:
#         values = {'hg_slope': hg_slope, 'hg_intercept': hg_intercept}
#         return(values[objName])
#     except(NameError, AttributeError):
#         return(np.NaN)



def MCI_calculation_full(acq, objName):
    '''Performs all the clumped isotope calculations for a single acq'''

    # calculating measured voltage ratios, with sample/ref bracketing
    R_measured_sample=(acq.voltSam[:,(4,7)]/np.tile(acq.voltSam[:,2],(2,1)).T)
    R_measured_ref=(acq.voltRef[:,(4,7)]/np.tile(acq.voltRef[:,2],(2,1)).T)

    # Subtracting adducts
    R_measured_sample[:,0] -= (acq.voltSam[:,2]*acq.adduct_17[1] + np.square(acq.voltSam[:,2])*acq.adduct_17[0])
    R_measured_sample[:,1] -= (acq.voltSam[:,4]*acq.adduct_18[1] + np.square(acq.voltSam[:,4])*acq.adduct_18[0])

    R_measured_ref[:,0] -= (acq.voltRef[:,2]*acq.adduct_17[1] + np.square(acq.voltRef[:,2])*acq.adduct_17[0])
    R_measured_ref[:,1] -= (acq.voltRef[:,4]*acq.adduct_18[1] + np.square(acq.voltRef[:,4])*acq.adduct_18[0])


    delta_measured=np.zeros(np.shape(R_measured_sample)) # Preallocating for size of delta array

    for l in range(len(R_measured_sample)):
        delta_measured[l,:]= (R_measured_sample[l,:]/((R_measured_ref[l,:]+R_measured_ref[l+1,:])/2)-1)*1000
        #correction for stretching
        delta_measured[l,1]= delta_measured[l,1]*np.mean(R_measured_ref[l:l+1,1])/acq.stretch_18

    # couches ratios in analysis/std bracketing, put in delta notation
    # correcting d18 measurement using fragmentation value
    delta_measured_mean=np.zeros((3,delta_measured.shape[1]))

    delta_measured_mean[0,:]=np.mean(delta_measured, axis=0) # averaging for all cycles
    delta_measured_mean[1,:]=np.std(delta_measured, axis=0,ddof=1)  # standard deviation among all cycles
    delta_measured_mean[2,:]=delta_measured_mean[1,:]/np.sqrt(len(delta_measured)) # std error among all cycles

    d17=delta_measured_mean[0,0]
    d18=delta_measured_mean[0,1]

    calculatedCIValues = {'d17_full': d17, 'd18': d18}

    return calculatedCIValues[objName]

def MCI_calculation_D(acq, objName):
    '''Performs all the clumped isotope calculations for a single acq'''

    # calculating measured voltage ratios, with sample/ref bracketing
    R_measured_sample = acq.voltSam[:,4]/acq.voltSam[:,2]
    R_measured_ref = acq.voltRef[:,4]/acq.voltRef[:,2]

    # Subtracting adducts
    R_measured_sample -= (acq.voltSam[:,2]*acq.adduct_17[1] + np.square(acq.voltSam[:,2])*acq.adduct_17[0])
    R_measured_ref -= (acq.voltRef[:,2]*acq.adduct_17[1] + np.square(acq.voltRef[:,2])*acq.adduct_17[0])

    delta_measured=np.zeros(np.shape(R_measured_sample)) # Preallocating for size of delta array

    for l in range(len(R_measured_sample)):
        delta_temp= (R_measured_sample[l]/((R_measured_ref[l]+R_measured_ref[l+1])/2)-1)*1000
        #correction for stretching
        delta_measured[l]= delta_temp*np.mean(R_measured_ref[l:l+1])/acq.stretch_17

    # couches ratios in analysis/std bracketing, put in delta notation
    # correcting d18 measurement using fragmentation value
    delta_measured_mean=np.zeros(3)

    delta_measured_mean[0]=np.mean(delta_measured, axis=0) # averaging for all cycles
    delta_measured_mean[1]=np.std(delta_measured, axis=0,ddof=1)  # standard deviation among all cycles
    delta_measured_mean[2]=delta_measured_mean[1]/np.sqrt(len(delta_measured)) # std error among all cycles

    d17=delta_measured_mean[0]

    calculatedCIValues = {'d17_D': d17}

    return calculatedCIValues[objName]

def MCI_bulk_solver(z,*extraArgs):
    '''Function for solving for bulk composition'''
    m, frag = extraArgs # unpacking tuple

    rD = 4*z[0]/(1+frag*(z[1]+3*z[0]))-m[0]
    r13C = (4*z[0]+z[1])/(1+frag*(z[1]+3*z[0]))-m[1]
    return(rD,r13C)

def MCI_r_to_c(rD,r13):
    ''' func for caculating methane stochastic isotopologue concentrations from ratios '''

    c12 = 1/(1+r13)
    c13 = r13/(1+r13)
    cH = 1/(1+rD)
    cD = rD/(1+rD)

    c12CH4 = c12*cH**4
    c12CH3D = 4*c12*cH**3*cD
    c13CH4 = c13*cH**4
    c13CH3D = 4*c13*cH**3*cD
    c12CH2D2 = 6*c12*cH**2*cD**2

    return(c12CH4, c12CH3D, c13CH4, c13CH3D, c12CH2D2)


def MCI_bulk_comp(analysis, objName):
    '''Performs all the clumped isotope calculations for a single acq'''

    r13_vpdb=0.0112372
    rD_vsmow=0.00015576

    r17_full = analysis.d17_full/1000 + 1
    d17_full_sterr = analysis.d17_full_sterr
    r17_D = analysis.d17_D/1000 + 1
    d17_D_sterr = analysis.d17_D_sterr
    r18 = analysis.d18/1000 + 1
    d18_sterr = analysis.d18_sterr

    rD_ref = (analysis.dD_ref/1000+1)*rD_vsmow
    r13_ref = (analysis.d13C_ref/1000+1)*r13_vpdb

    c12CH4_ref, c12CH3D_ref, c13CH4_ref, c13CH3D_ref, c12CH2D2_ref = MCI_r_to_c(rD_ref,r13_ref)

    m_1 = r17_D*(4*rD_ref)/(1+analysis.frag*(r13_ref+3*rD_ref))
    m_2 = r17_full*(4*rD_ref+r13_ref)/(1+analysis.frag*(r13_ref+3*rD_ref))
    m_3 = r18*(6*rD_ref**2+4*rD_ref*r13_ref)/(1+analysis.frag*(r13_ref+3*rD_ref))

    m = [m_1, m_2]
    z_guess = [0, 0]
    extraArgs = (m, analysis.frag,)

    bulkSolverResult = root(MCI_bulk_solver,z_guess, args = extraArgs, method = 'lm', options = {'xtol':1e-20, 'ftol':1e-20})

    rD_sa, r13_sa = bulkSolverResult.x

    c12CH4_sa, c12CH3D_sa, c13CH4_sa, c13CH3D_sa, c12CH2D2_sa = MCI_r_to_c(rD_sa,r13_sa)

    D18 = 1000*(m_3*(c12CH4_sa + analysis.frag*(0.75*c12CH3D_sa+c13CH4_sa))/(c13CH3D_sa+c12CH2D2_sa)-1)

    dD_sa = (rD_sa/rD_vsmow-1)*1000
    d13_sa = (r13_sa/r13_vpdb-1)*1000


    calculatedMCIValues = {'dD': dD_sa, 'd13C': d13_sa, 'D18_raw': D18, 'dD_sterr': analysis.d17_D_sterr, 'd13C_sterr': analysis.d17_full_sterr, 'D18_sterr': analysis.d18_sterr}

    return calculatedMCIValues[objName]

def MCI_averages(analysis, objName):
    '''Computes the mean, std dev, and std error for every attribute useful for a CI measurement'''

    props2=['d17_full','d17_full_sterr','d17_D','d17_D_sterr','d18','d18_sterr']

    valName = objName.replace('_sterr','')
    values=[]#preallocating for value storage

    if valName == 'd17_D':
        acqsToUse = range(len(analysis.acqs_D))
        for i in acqsToUse:
            values.append(getattr(analysis.acqs_D[i],valName))
    else:
        acqsToUse = range(len(analysis.acqs_full))
        for i in acqsToUse:
            values.append(getattr(analysis.acqs_full[i],valName))

    if len(values) != 0:
        values = np.asarray(values)
        if '_stdev' in objName:
            return values.std(axis=0,ddof=1)
        elif '_sterr' in objName:
            return values.std(axis=0,ddof=1)/np.sqrt(values.shape[0])
        else:
            return values.mean(axis=0)


def Sample_type_checker(analyses):
    AllAnalysesHaveType = True
    for i in range(len(analyses)):
        if analyses[i].type not in ['std1','sample','stdD13','hg']:
            AllAnalysesHaveType = False
            break

    return AllAnalysesHaveType

def Get_types_auto(analyses):
    '''Function to assign the correct type to every analysis automatically, given a few assumptions:
    stds are +1 std or +D+13, heated gases have ' hg ' in name, everythin else is a sample)'''

    print('Automatically assigning analyses types ')
    choice = raw_input('(s)top process, see (n)aming guidelines, or hit any other key to continue ').lower()
    if choice == 's':
        return(analyses)
    elif choice == 'n':
        print('1. +1 std must contain "std" AND "+1" ')
        print('2. +D+13 std "+D" AND "+13" ')
        print('3. Heated gases must contain "hg" AND < 14 chars ')
        print('4. Analyses that already have a valid type are not modified ')
        print('')
        return(analyses)
    else:
        for i in range(len(analyses)):
            if analyses[i].type in ['std1', 'stdD13', 'hg', 'sample']:
                continue
            else:
                name = analyses[i].name.lower()
                if ('+1' in name) and ('std' in name) and ('+d' not in name):
                    analyses[i].name = str('"' + analyses[i].name + '"')
                    analyses[i].type = 'std1';
                    continue
                elif ('+d+13' in name):
                    analyses[i].type = 'stdD13'
                    analyses[i].name = str('"' + analyses[i].name + '"')
                    continue
                elif ('hg' in name) and (len(name) < 14):
                    analyses[i].type = 'hg'
                    continue
                else:
                    analyses[i].type = 'sample'
    if Sample_type_checker(analyses):
        print('Types successfully assigned ')
    else:
        print('Automatic assignment failed ')

    return(analyses)

def Get_types_manual(analyses):
    '''Function to assign the correct type to every analysis manually'''

    print('Manually assigning analyses types ')
    for i in range(len(analyses)):
        if analyses[i].type not in ['std1', 'stdD13', 'sample', 'hg']:
            typeChoice = raw_input('Type for: ' + analyses[i].name + ' -> (h)g, (s)ample, std(D)13, or std(1)?').lower()
            if typeChoice == 'd':
                analyses[i].type = 'stdD13'
                analyses[i].name = str('"' + analyses[i].name + '"')
            elif typeChoice == 'h':
                analyses[i].type = 'hg'
            elif typeChoice == '1':
                analyses[i].type = 'std1'
                analyses[i].name = str('"' + analyses[i].name + '"')
            else:
                analyses[i].type = 'sample'

    if Sample_type_checker(analyses):
        print('Types successfully assigned ')
    else:
        print('Assignment of some types failed ')

    return(analyses)


def ExportSequence(analyses):
    '''Most common export sequence'''
    print('Exporting to temporary PAMS and FlatList files ')
    print('Exporting full acqs to a PAMS sheet...')
    exportNamePAMS = 'autoPAMS_Export'
    PAMS_exporter(analyses, exportNamePAMS, displayProgress = True)
    print('Exporting analyses to a flatlist...')
    exportNameFlatlist = 'autoFlatListExport'
    FlatList_exporter(analyses,exportNameFlatlist, displayProgress = True)
    print('Analyses successfully exported')

    return


def MCI_background_correction(instance, objName):
    # voltSamTemp = np.copy(instance.voltSam_raw)
    # voltRefTemp = np.copy(instance.voltRef_raw)

    backgroundArray = np.array([0, 0 , instance.background[0], 0, instance.background[1], 0, 0, instance.background[2], 0, 0])

    # absolute background correction
    try:


        voltSamTemp = instance.voltSam_raw - np.tile(backgroundArray,(len(instance.voltSam_raw),1))
        voltRefTemp = instance.voltRef_raw - np.tile(backgroundArray,(len(instance.voltRef_raw),1))
    except(IndexError):
        return 0

    voltages = {'voltSam': voltSamTemp, 'voltRef': voltRefTemp}

    return(voltages[objName])


def lsqfitma(X, Y):
    """
        From:
        # teaching.py
    #
    # purpose:  Teaching module of ff_tools.
    # author:   Filipe P. A. Fernandes
    # e-mail:   ocefpaf@gmail
    # web:      http://ocefpaf.tiddlyspot.com/
    # created:  09-Sep-2011
    # modified: Sun 23 Jun 2013 04:30:45 PM BRT
    #
    Calculate a "MODEL-2" least squares fit.

    The line is fit by MINIMIZING the NORMAL deviates.

    The equation of the line is:     y = mx + b.

    This line is called the MAJOR AXIS.  All points are given EQUAL
    weight.  The units and range for X and Y must be the same.
    Equations are from York (1966) Canad. J. Phys. 44: 1079-1086;
    re-written from Kermack & Haldane (1950) Biometrika 37: 30-41;
    after a derivation by Pearson (1901) Phil. Mag. V2(6): 559-572.

    Data are input and output as follows:

    m, b, r, sm, sb = lsqfitma(X, Y)
    X    =    x data (vector)
    Y    =    y data (vector)
    m    =    slope
    b    =    y-intercept
    r    =    correlation coefficient
    sm   =    standard deviation of the slope
    sb   =    standard deviation of the y-intercept

    Note that the equation passes through the centroid:  (x-mean, y-mean)

    """

    X, Y = map(np.asanyarray, (X, Y))

    # Determine the size of the vector.
    n = len(X)

    # Calculate sums and other re-used expressions.
    Sx = np.sum(X)
    Sy = np.sum(Y)
    xbar = Sx / n
    ybar = Sy / n
    U = X - xbar
    V = Y - ybar

    Suv = np.sum(U * V)
    Su2 = np.sum(U ** 2)
    Sv2 = np.sum(V ** 2)

    sigx = np.sqrt(Su2 / (n - 1))
    sigy = np.sqrt(Sv2 / (n - 1))

    # Calculate m, b, r, sm, and sb.
    m = (Sv2 - Su2 + np.sqrt(((Sv2 - Su2) ** 2) + (4 * Suv ** 2))) / (2 * Suv)
    b = ybar - m * xbar
    r = Suv / np.sqrt(Su2 * Sv2)

    sm = (m / r) * np.sqrt((1 - r ** 2) / n)
    sb1 = (sigy - sigx * m) ** 2
    sb2 = (2 * sigx * sigy) + ((xbar ** 2 * m * (1 + r)) / r ** 2)
    sb = np.sqrt((sb1 + ((1 - r) * m * sb2)) / n)

    return m, b, r, sm, sb


def lsqcubic(X, Y, sX, sY, tl=1e-7):
    """
    From:
    # teaching.py
#
# purpose:  Teaching module of ff_tools.
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  09-Sep-2011
# modified: Sun 23 Jun 2013 04:30:45 PM BRT
#
# obs: Just some basic example function.

    Calculate a MODEL-2 least squares fit from weighted data.

    The line is fit by MINIMIZING the weighted residuals in both x & y.
    The equation of the line is:     y = mx + b,
    where m is determined by finding the roots to the cubic equation:

    m^3 + P * m^2 + Q * m + R = 0.

    Eqs for P, Q and R are from York (1966) Canad. J. Phys. 44: 1079-1086.

    Data are input and output as follows:
    m, b, r, sm, sb, xc, yc, ct = lsqcubic(X, Y, sX, sY, tl)
    X    =    x data (vector)
    Y    =    y data (vector)
    sX   =    uncertainty of x data (vector)
    sY   =    uncertainty of y data (vector)
    tl   =    test limit for difference between slope iterations

    m    =    slope
    b    =    y-intercept
    r    =    weighted correlation coefficient
    sm   =    standard deviation of the slope
    sb   =    standard deviation of the y-intercept
    xc   =    WEIGHTED mean of x values
    yc   =    WEIGHTED mean of y values
    ct   =    count: number of iterations

    Notes:  1.  (xc,yc) is the WEIGHTED centroid.
            2.  Iteration of slope continues until successive differences
                are less than the user-set limit "tl".  Smaller values of
                tl require more iterations to find the slope.
            3.  Suggested values of tl = 1e-4 to 1e-6.

    """

    X, Y = map(np.asanyarray, (X, Y))

    # Find the number of data points and make one time calculations:
    n = len(X)
    wX = 1 / (sX ** 2)
    wY = 1 / (sY ** 2)

    # Set-up a few initial conditions:
    ct, ML = 0, 1

    # ESTIMATE the slope by calculating the major axis according
    # to Pearson's (1901) derivation, see: lsqfitma.

    MC = lsqfitma(X, Y)[0]

    test = np.abs((ML - MC) / ML)

    # Calculate the least-squares-cubic. Make iterative calculations until the
    # relative difference is less than the test conditions

    while test > tl:
        # Calculate sums and other re-used expressions:
        MC2 = MC ** 2
        W = (wX * wY) / ((MC2 * wY) + wX)
        W2 = W ** 2

        SW = np.sum(W)
        xc = (np.sum(W * X)) / SW
        yc = (np.sum(W * Y)) / SW

        U = X - xc
        V = Y - yc

        U2 = U ** 2
        V2 = V ** 2

        SW2U2wX = np.sum(W2 * U2 / wX)

        # Calculate coefficients for least-squares cubic:
        P = -2 * np.sum(W2 * U * V / wX) / SW2U2wX
        Q = (np.sum(W2 * V2 / wX) - np.sum(W * U2)) / SW2U2wX
        R = np.sum(W * U * V) / SW2U2wX
        # Find the roots to the least-squares cubic:
        LSC = [1, P, Q, R]
        MR = np.roots(LSC)

        # Find the root closest to the slope:
        DIF = np.abs(MR - MC)
        MinDif, Index = DIF.min(), DIF.argmin()

        ML = MC
        MC = MR[Index]
        test = np.abs((ML - MC) / ML)
        ct = ct + 1

    # Calculate m, b, r, sm, and sb.
    m = MC
    b = yc - m * xc
    r = np.sum(U * V) / np.sqrt(np.sum(U2) * np.sum(V2))
    sm2 = (1 / (n - 2)) * (np.sum(W * (((m * U) - V) ** 2)) / np.sum(W * U2))
    sm = np.sqrt(sm2)
    sb = np.sqrt(sm2 * (np.sum(W * (X ** 2)) / SW))

    return m, b, r, sm, sb, xc, yc, ct
