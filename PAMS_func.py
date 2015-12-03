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
from scipy.optimize import root


class MCI_VALUE(object):
    '''subclass defining how important isotopic ratios are calculated, and d18O_mineral'''

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
        if len(instance.acqs)>=1:
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
        if self.name in ['d13C', 'd13C_sterr', 'dD', 'dD_sterr', 'D18', 'D18_sterr']:
            return np.around(MCI_bulk_comp(instance, self.name),3)

        else:
            raise ValueError('Sample D47_raw is out of range')


    def __set__(self, obj, value):
        raise AttributeError('Cannot change CI correction scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete CI corretion scheme')

# class CI_UNIVERSAL_VALUE(object):
#     '''subclass defining how heated gases are taken'''
#
#     def __init__(self, name):
#         self.name = name
#
#     def __get__(self,instance,cls):
#         if self.name in ['hg_slope', 'hg_intercept']:
#             return np.around(CI_hg_values(instance, self.name),5)
#
#         else:
#             raise ValueError('Not a valid value')
#
#     def __set__(self, obj, value):
#         raise AttributeError('Cannot change individual hg value')
#
#     def __delete__(self, instance):
#         raise AttributeError('Cannot delete individual hg value')

class MCI_VOLTAGE(object):
    '''subclass defining how raw voltages are corrected'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if self.name in ['voltSam', 'voltRef']:
            return np.around(CI_background_correction(instance, self.name),3)

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

        self.full_adduct = [0, 0]
        self.D_adduct = [0, 0]
        self.adduct_18 = [0, 0]

        self.mass18_frag = 3.700725
        self.mass17_frag = 0.053245
        self.frag = 0.789



    d17_full = MCI_AVERAGE('d17_full')
    d18 = MCI_AVERAGE('d18')
    d17_D = MCI_AVERAGE('d17_D')

    D18_raw = MCI_CALCULATED_VALUE('D18_raw')
    d13C = MCI_CALCULATED_VALUE('d13C')
    dD = MCI_CALCULATED_VALUE('dD')

    d17_full_sterr = MCI_AVERAGE('d17_full_sterr')
    d18_sterr = MCI_AVERAGE('d18_sterr')
    d17_D_sterr = MCI_AVERAGE('d17_full_sterr')

    D18_sterr = MCI_CALCULATED_VALUE('D18_sterr')
    d13C_sterr = MCI_CALCULATED_VALUE('d13C_sterr')
    dD_sterr = MCI_CALCULATED_VALUE('dD_sterr')

    # D47_CRF = CI_CORRECTED_VALUE('D47_CRF')
    # hg_slope = CI_UNIVERSAL_VALUE('hg_slope')
    # hg_intercept = CI_UNIVERSAL_VALUE('hg_intercept')



class ACQUISITION_FULL(object):
    "A class for all the attributes of a single clumped isotope acquision"


    def __init__(self,acqNum):
        self.acqNum=acqNum
        self.voltSam_raw=[]
        self.voltRef_raw=[]
        self.background=[]
        self.date=''
        self.time=''
        self.adduct_17 = [ 0, 0]
        self.adduct_18 = [ 0, 0]
        self.frag_18 = 3.700725
        self.type = 'full'
        self.name

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
        self.background=[]
        self.date=''
        self.time=''
        self.adduct_17 = [0, 0]
        self.frag_17= 0.053245
        self.type = 'D'

    d17_D=MCI_VALUE('d17_D')
    voltSam = MCI_VOLTAGE('voltSam')
    voltRef = MCI_VOLTAGE('voltRef')



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
    try:
        analysisName = nameBlock[(nameBlock.find('Pressadjust')+19):(nameBlock.find('Identifier')-2)]
        # analysisName = filter(lambda)
    except UnicodeEncodeError:
        print('caught a unicode encode error')
    #quick function to filter out non-ascii chars that are tacked on to end
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



    return voltRef_raw, voltSam_raw, analysisName, date_str, time_str




def CIDS_cleaner(analyses):
    '''function for cleaning up a parsed CIDS file and alerting to any corrupted analyses'''

    # Cleaning up the parsed file
    # because of the way the parser is written, an extra empty analysis is added to the end of the list
    if not analyses[-1].acqs:
        del analyses[-1]


    # checking for analyses with not enough acqs, alerting if there are some
    temp=8
    lowAcqs= [k for k in analyses if len(k.acqs)<temp]
    if not lowAcqs:
        print 'All analyses have at least %d acquistions' % temp

    else:
        print 'analyses with too few acqs are:'
        print '\n'.join(['Sample '+ str(k.num)+ ': '+ k.name +' has ' + str(len(k.acqs)) + ' acquisitions' for k in lowAcqs])

    # converting the voltages and backgrounds to arrays so that they can be easily used to do calculations
    for i in range(len(analyses)):
            for j in range(len(analyses[i].acqs)):
                # converting voltagted to arrays, and doing 3 sig figs to match CIDS Sheet
                analyses[i].acqs[j].voltSam_raw=np.around(np.asarray(analyses[i].acqs[j].voltSam_raw),3)
                analyses[i].acqs[j].voltRef_raw=np.around(np.asarray(analyses[i].acqs[j].voltRef_raw),3)
                analyses[i].acqs[j].background=np.asarray(analyses[i].acqs[j].background)
                # rounding d13C and d18O of each acq to 3 sig figs to match CIDS sheet
                (analyses[i].acqs[j].d13C,analyses[i].acqs[j].d18O_gas) = np.around((analyses[i].acqs[j].d13C,analyses[i].acqs[j].d18O_gas),3)


    print 'All analyses are cleaned, and voltages converted to arrays'

    return analyses

def FlatList_exporter(analyses,fileName, displayProgress = False):
    '''Exports a CSV file with the essential calculated MCI info'''

    export=open(fileName + '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    wrt.writerow(['User','date','Type','Sample ID','spec #\'s', 'full acqs', 'dD acqs', 'd13C (vpdb)','d13C_sterr','d2H (vsmow)',
    'd2H_sterr','d18','d18_sterr','D18 (v. wg)','D18_sterr'])
    counter = 0
    if displayProgress:
        for item in analyses:
            wrt.writerow([item.user, item.date, item.type, item.name, item.num, len(item.acqs_full), len(item.acqs_D), item.d13C, item.d13C_stdev, item.dD,
            item.dD_sterr,item.d18,item.d18_sterr,item.D18,item.D18_sterr ])
            counter += 1
            if ((counter * 100)*100) % (len(analyses)*100) == 0:
                print(str((counter*100)/len(analyses)) + '% done')
    else:
        for item in analyses:
            wrt.writerow([item.user, item.date, item.type, item.name, item.num, len(item.acqs_full), len(item.acqs_D), item.d13C, item.d13C_stdev, item.dD,
            item.dD_sterr,item.d18,item.d18_sterr,item.D18,item.D18_sterr])
    export.close()
    return

def CIDS_exporter(analyses, fileName, displayProgress = False):
    '''Exports a CSV file that is the same format as a traditional CIDS file.
    Python CIDS files are also importable, to load all important sample info
    back into the program'''

    export=open(fileName + '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    for analysis in analyses:
        wrt.writerow(['__NewSample__'])
        for acq in analysis.acqs_D:
            wrt.writerow(['__NewAcq-D__'])
            wrt.writerow(['__AcqVoltage__','',]+np.shape(acq.voltSam)[1]*['',]+acq.voltRef[0,:].tolist())
            for i in range(len(acq.voltSam)):
                wrt.writerow(['',''] + acq.voltSam[i,:].tolist() + acq.voltRef[i+1,:].tolist())
            wrt.writerow([])
            wrt.writerow(['','acq num','date','name','type', 'adductLine (2nd-ordr)', 'adducLine (1st-ordr)','17fragCorrection','d17_D (wg)'])
            wrt.writerow(['__AcqInfo__',acq.acqNum, acq.date, analysis.name, acq.type, acq.adduct_17[0], acq.adduct_17[1], acq.frag_17, acq.d17_D])

        wrt.writerow([])

        for acq in analysis.acqs_full:
            wrt.writerow(['__NewAcq-Full__'])
            wrt.writerow(['__AcqVoltage__','',]+np.shape(acq.voltSam)[1]*['',]+acq.voltRef[0,:].tolist())
            for i in range(len(acq.voltSam)):
                wrt.writerow(['',''] + acq.voltSam[i,:].tolist() + acq.voltRef[i+1,:].tolist())
            wrt.writerow([])
            wrt.writerow(['','acq num','date','name','type', '17adductLine (2nd-ordr)',
            '17adductLine (1st-ordr)', '18adductLine (2nd-ordr)', '18adductLine (1st-ordr)','18fragCorrection','d17_full (wg)', 'd18 (wg)'])
            wrt.writerow(['__AcqInfo__',acq.acqNum, acq.date, analysis.name, acq.type, acq.adduct_17[0], acq.adduct_17[1],
            acq.adduct_18[0], acq.adduct_18[1],acq.frag_18, acq.d17_full, acq.d18])

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
        wrt.writerow(['','User','date','Type','Sample ID','acq #\'s','D_adduct (2nd-ordr)',
        'D_adduct (1st-ordr)','full_adduct (2nd-ordr)','full_adduct (1st-ordr)','18_adduct (2nd-ordr)',
        '18_adduct (1st-ordr)','D_fragCorrection','18_fragCorrection','Fragmenation','d13C (vpdb)',
        'd13C_sterr','dD (vsmow)','dD_sterr','d47','d47_stdev','D47 (v. Oz)','D47_stdev','D47_sterr','d48', 'd48_stdev','D48','D48_stdev'])
        wrt.writerow(['__SampleSummary__',analysis.user, analysis.date, analysis.type, analysis.name, analysis.num, analysis.skipFirstAcq, analysis.d13C, analysis.d13C_stdev,
        analysis.d18O_gas, analysis.d18O_min, analysis.d18O_stdev, analysis.d47, analysis.d47_stdev, analysis.D47_raw,
        analysis.D47_stdev, analysis.D47_sterr, analysis.d48, analysis.d48_stdev, analysis.D48_raw, analysis.D48_stdev])
        wrt.writerow(24*['---',])

    export.close()
    return

def CIDS_importer(filePath, displayProgress = False):
    '''Imports voltage and analysis data that were exported using the CIDS_exporter function'''
    fileImport = open(filePath,'rU')
    fileReader = csv.reader(fileImport, dialect='excel')
    analyses = []

    cycles_in_acqs = 7

    for line in fileReader:
        if '__NewSample__' in line:
            analyses.append(CI())
            continue
        if '__NewAcq__' in line:
            analyses[-1].acqs.append(ACQUISITION(0))
            continue
        if '__AcqVoltage__' in line:
            voltIndex = line.index('__AcqVoltage__')
            cycle = 0
            voltRefs = [float(vr) for vr in line[voltIndex+8:voltIndex+14]]
            analyses[-1].acqs[-1].voltRef_raw.append(voltRefs)
            continue
        if cycle < (cycles_in_acqs):
            voltSams=[float(vs) for vs in line[voltIndex+2:voltIndex+8]]
            voltRefs = [float(vr) for vr in line[voltIndex+8:voltIndex+14]]
            analyses[-1].acqs[-1].voltSam_raw.append(voltSams)
            analyses[-1].acqs[-1].voltRef_raw.append(voltRefs)
            cycle += 1
            continue
        if '__AcqInfo__' in line:
            acqIndex = line.index('__AcqInfo__')
            analyses[-1].acqs[-1].acqNum = float(line[acqIndex + 1])
            analyses[-1].acqs[-1].date = line[acqIndex + 2]
            analyses[-1].acqs[-1].d13Cref = float(line[acqIndex + 3])
            analyses[-1].acqs[-1].d18Oref = float(line[acqIndex + 4])
            analyses[-1].acqs[-1].d13C = float(line[acqIndex + 7])
            analyses[-1].acqs[-1].d18O_gas = float(line[acqIndex + 8])
            analyses[-1].acqs[-1].voltRef_raw = np.asarray(analyses[-1].acqs[-1].voltRef_raw)
            analyses[-1].acqs[-1].voltSam_raw = np.asarray(analyses[-1].acqs[-1].voltSam_raw)
            continue
        if '__SampleSummary__' in line:
            summaryIndex = line.index('__SampleSummary__')
            [analyses[-1].user, analyses[-1].date, analyses[-1].type, analyses[-1].name] = line[summaryIndex + 1:summaryIndex+5]
            analyses[-1].num = int(line[summaryIndex + 5])
            analyses[-1].skipFirstAcq = bool(line[summaryIndex + 6])

    fileImport.close()
    return analyses




def Get_gases(analyses):
    '''Finds which analyses are heated and equilibrated gases, and assigns them TCO2 values'''
    properNames = raw_input("Do all equilibrated gases have '25' in name? (y/n) ").lower()
    for item in analyses:
        if 'BOC' in item.name.upper():
            if properNames == 'y':
                if '25' in item.name:
                    item.TCO2 = 25
                    item.type = 'eg'
                else :
                    item.TCO2 = 1000
                    item.type = 'hg'
            else :
                while True:
                    choice = raw_input('Is acq num: ' + str(item.num) + ' with name: ' + item.name + ' a (h)eated gas, an (e)quilibrated gas?, or (s)kip? ')
                    if choice.lower() == 'h':
                        item.TCO2 = 1000
                        item.type = 'hg'
                        break
                    elif choice.lower() == 'e':
                        item.TCO2 = 25
                        item.type = 'eg'
                        break
                    elif choice.lower() == 's':
                        print('Skipping this analysis ')
                        break
                    else:
                        print('Not a valid entry, please try again')
    print('All heated and equilibrated gases have been found and labeled ')
    return

def Daeron_exporter(analyses, fileName):
    '''Exports analyses in a csv that is formatted for Matthieu Daeron's xlcor47 script'''

    export=open(fileName +'_daeron'+ '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    for item in analyses:
        if np.isnan(item.TCO2):
            wrt.writerow([item.name, item.d47, item.D47_raw, item.D47_sterr, ])
        else :
            wrt.writerow([item.name, item.d47, item.D47_raw, item.D47_sterr, item.TCO2])
    export.close()
    return


def Pressure_Baseline_Processer(fileFolder):
    '''Pulls raw intensity data out of a folder of peak scans'''

    fileList = [i for i in os.listdir(fileFolder) if '.csv' in i] #copies all csv pbl files in directory
    print fileList
    day = '4-22-14' # TODO: make this I/O eventually
    fileList = [i for i in fileList if day in i] #only the csv files from certain days
    A=[]
    print fileList

    #First read in the csv files to numpy arrays
    for i in fileList:
        A.append([])
        fileName = fileFolder+i
        fileReader = csv.reader(open(fileName, 'rU'), dialect='excel')
        for line in fileReader:
            A[-1].append(line)
        A[-1].pop(0) #remove the first line with the headers

        A[-1]=np.asfarray(A[-1])
        # as a note, columns, in order, are: scanNumber, HV, mass44, mass45, mass46, mass47, mass48, mass49

    # Now, finding minimums
    int44=np.zeros((len(A),1))
    minimums=np.zeros((len(A),6))
    int49=np.zeros((len(A),1))

    for i in range(len(A)):
        max44=A[i][:,2].max()
        halfHeight=np.nonzero(A[i][:,2]>(max44/2)) #array of indices of all values with mass44 intensity is gt half of max
        peakCenter=int(np.median(halfHeight))
        halfWidth=halfHeight[0][-1]-peakCenter # currently this line is not necessary, but may be later on
        int44[i] = np.mean(A[i][peakCenter-5:peakCenter+5,2])
        int49[i] = np.mean(A[i][peakCenter-5:peakCenter+5,7])
        minimums[i,:]=A[i][peakCenter:,2:].min(axis=0) # minimum values of mass 44-49 on the right side of the peak





    return int44, int49, minimums

def MCI_calculation_full(acq, objName):
    '''Performs all the clumped isotope calculations for a single acq'''

    # calculating measured voltage ratios, with sample/ref bracketing
    R_measured_sample=(acq.voltSam[:,(4,7)]/np.tile(acq.voltSam[:,2],(2,1)).T)
    R_measured_ref=(acq.voltRef[:,(4,7)]/np.tile(acq.voltRef[:,2],(2,1)).T)

    # Subtracting adducts
    R_measured_sample[:,0] -= (voltSam[:,2]*adduct_17[0] + np.square(voltSam[:,2])*adduct_17[1])
    R_measured_sample[:,1] -= (voltSam[:,2]*adduct_18[0] + np.square(voltSam[:,2])*adduct_18[1])

    R_measured_ref[:,0] -= (voltRef[:,2]*adduct_17[0] + np.square(voltRef[:,2])*adduct_17[1])
    R_measured_ref[:,1] -= (voltRef[:,2]*adduct_18[0] + np.square(voltRef[:,2])*adduct_18[1])


    delta_measured=np.zeros(np.shape(R_measured_sample)) # Preallocating for size of delta array

    for l in range(len(R_measured_sample)):
        delta_measured[l,:]= (R_measured_sample[l,:]/((R_measured_ref[l,:]+R_measured_ref[l+1,:])/2)-1)*1000
        #correction for fragmentation
        delta_measured[l,1]= delta_measured[l,1]*np.mean(R_measured_ref[l:l+1,1])/acq.frag_18

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
    R_measured_sample -= (voltSam[:,2]*adduct_17[0] + np.square(voltSam[:,2])*adduct_17[1])
    R_measured_ref -= (voltRef[:,2]*adduct_17[0] + np.square(voltRef[:,2])*adduct_17[1])

    delta_measured=np.zeros(np.shape(R_measured_sample)) # Preallocating for size of delta array

    for l in range(len(R_measured_sample)):
        delta_temp= (R_measured_sample[l]/((R_measured_ref[l]+R_measured_ref[l+1])/2)-1)*1000
        #correction for fragmentation
        delta_measured[l]= delta_temp*np.mean(R_measured_ref[l:l+1])/acq.frag_17

    # couches ratios in analysis/std bracketing, put in delta notation
    # correcting d18 measurement using fragmentation value
    delta_measured_mean=np.zeros(3)

    delta_measured_mean[0]=np.mean(delta_measured, axis=0) # averaging for all cycles
    delta_measured_mean[1]=np.std(delta_measured, axis=0,ddof=1)  # standard deviation among all cycles
    delta_measured_mean[2]=delta_measured_mean[1]/np.sqrt(len(delta_measured)) # std error among all cycles

    d17=delta_measured_mean[0]

    calculatedCIValues = {'d17_full': d17}

    return calculatedCIValues[objName]

def MCI_bulk_solver(z,*extraArgs):
    '''Function for solving for bulk composition'''
    m, frag = extraArgs # unpacking tuple

    rD = 4*z[0]/(1+frag*(z[1]+3*z[0]))-m[0]
    r13C = (4*z[0]+z[1])/(1+frag*(z[1]+3*z[0]))-m[0]
    return(rD,r13C)

def MCI_r_to_c(rD,r13):
    ''' func for caculating methane stochastic isotopologue concentrations from ratios '''

    c12 = 1/(1+r13)
    c13 = r13/(1+r13)
    cH = 1/(1+rD)
    cD = rD/(1+rD)

    c12CH4 = c12*cH^4
    c12CH3D = 4*c12*cH^3*cD
    c13CH4 = c13*cH^4
    c13CH3D = 4*c13*cH^3*cD
    c12CH2D2 = 6*c12*cH^2*cD^2

    return(c12CH4, c12CH3D, c13CH4, c13CH3D, C12CH2D2)


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

    m_1 = r17_D*(4*rD_ref)/(1+analysis.frag*(r13C+3*rD_ref))
    m_2 = r17_full*(4*rD_ref+r13_ref)/(1+analysis.frag*(r13_ref+3*rD_ref))
    m_3 = r18*(6*rD_ref^2+4*rD_ref*r13_ref)/(1+analysis.frag*(r13_ref+3*rD_ref))

    m = [m_1, m_2]
    z_guess = [0, 0]
    extraArgs = (m, analysis.frag)

    [rD_sa, r13_sa] = root(MCI_bulk_solver,z_guess, args = extraArgs, method = 'lm', options = {'xtol':1e-20, 'ftol':1e-20})

    c12CH4_sa, c12CH3D_sa, c13CH4_sa, c13CH3D_sa, c12CH2D2_sa = MCI_r_to_c(rD_sa,r13_sa)

    D18 = 1000*(m_3*c12CH4 + analysis.frag(0.75*c12CH3D_sa+c13CH4_sa)/(c13CH3D_sa+c12CH2D2_sa)-1)

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

def CI_comparer(analyses1, analyses2):
    '''compares two sets of CI data for equivalence in D47 calculations'''
    equalSoFar = True
    if not len(analyses1) == len(analyses2):
        print('Data sets are not the same length!')

    for i in range(len(analyses1)):
        for j in range(len(analyses1[i].acqs)):
            if not analyses1[i].acqs[j].D47_raw == analyses2[i].acqs[j].D47_raw:
                print('Acq num {1} from analysis number {0} from comparison does not agree in D47'.format(i,j))
                print('analyses1 D47 = {0:.3f}, analyses2 D47 = {1:.3f}'.format(analyses1[i].acqs[j].D47_raw, analyses2[i].acqs[j].D47_raw))
                equalSoFar = False
        if not analyses1[i].D47_stdev == analyses2[i].D47_stdev:
            print('Sample number {0} from data sets do not agree in D47_stdev'.format(i))
            print('1. D47_stdev = {0:.3f}, 2. D47_stdev = {0:.3f}'.format(analyses1[i].D47_stdev, analyses2[i].D48_stdev))

    if equalSoFar:
        print('Data sets are equivalent')

def Sample_type_checker(analyses):
    AllAnalysesHaveType = True
    for i in range(len(analyses)):
        if analyses[i].type not in ['std','sample','eg','hg']:
            AllAnalysesHaveType = False
            break

    return AllAnalysesHaveType

def Get_types_auto(analyses):
    '''Function to assign the correct type to every analysis automatically, given a few assumptions:
    stds are Carrara, NBS-19, or TV03, all gases have 'BOC' in name, and all egs have '25' in name)'''

    print('Automatically assigning analyses types ')
    choice = raw_input('(s)top process, see (n)aming guidelines, or hit any other key to continue ').lower()
    if choice == 's':
        return(analyses)
    elif choice == 'n':
        print('1. Standards must contain words "carrara", "TV03", or "NBS-19" ')
        print('2. 25 C equilibrated gases must contain "BOC" AND "25" ')
        print('3. 1000 C heated gases must contain "BOC" AND NOT "25" AND < 10 chars ')
        print('4. Analyses that already have a valid type are not modified ')
        return(analyses)
    else:
        for i in range(len(analyses)):
            if analyses[i].type in ['std', 'eg', 'hg', 'sample']:
                continue
            else:
                name = analyses[i].name.lower()
                if ('carrara' in name) or ('tv03' in name) or ('nbs-19' in name):
                    analyses[i].type = 'std';
                    continue
                elif ('boc' in name):
                    if ('25' in name):
                        analyses[i].type = 'eg'
                        continue
                    elif (len(name) < 10):
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
        if analyses[i].type not in ['eg', 'hg', 'sample', 'std']:
            typeChoice = raw_input('Type for: ' + analyses[i].name + ' -> (e)g, (h)g, (s)ample, or s(t)d?').lower()
            if typeChoice == 'e':
                analyses[i].type = 'eg'
            elif typeChoice == 'h':
                analyses[i].type = 'hg'
            elif typeChoice == 't':
                analyses[i].type = 'std'
            else:
                analyses[i].type = 'sample'

    if Sample_type_checker(analyses):
        print('Types successfully assigned ')
    else:
        print('Assignment of some types failed ')

    return(analyses)

def CI_CRF_data_corrector(analyses, showFigures = True):
    '''Extended function to do all aspects of CRF data correction'''

    # First, make a list of all hgs

    # Next, make a hg D48 line, using a York regression
    # hg_slope_48, hg_intercept_48, r_48, sm, sb = lsqfitma([i.d48 for i in hgs], [i.D48_raw for i in hgs])
    CI_48_excess_checker(analyses)
    # Make a new hg collection with any 48 excesses removed
    hgs = [i for i in analyses if (i.type == 'hg' and not i.D48_excess)]
    d47_hgs = np.asarray([i.d47 for i in hgs])
    D47_raw_hgs = np.asarray([i.D47_raw for i in hgs])
    d47_stdev_hgs = np.asarray([i.d47_stdev for i in hgs])
    D47_sterr_hgs = np.asarray([i.D47_sterr for i in hgs])
    # Now, make hg D47 line, using a York regression
    global hg_slope, hg_intercept
    hg_slope, hg_intercept, r_47, sm, sb, xc, yc, ct = lsqcubic(d47_hgs, D47_raw_hgs,d47_stdev_hgs, D47_sterr_hgs)

    D47_raw_hgs_model = d47_hgs*hg_slope+hg_intercept
    if showFigures:
        plt.figure(0)
        plt.figure(0).hold(True)
        plt.errorbar(d47_hgs, D47_raw_hgs, xerr = d47_stdev_hgs, yerr = D47_sterr_hgs, fmt = 'bo')
        plt.plot(d47_hgs, D47_raw_hgs_model, '-')
        plt.ylabel(ur'$\Delta_{47} \/ (\u2030)$')
        plt.xlabel(ur'$\delta^{47} \/ (\u2030)$')
        plt.show()

        # Calculation of the D47_CRF values are done automatically given its definition, based on the global hg slope and int
        # Plot the stds
        stds_Carrara = [i for i in analyses if (i.type == 'std' and 'carrara' in i.name.lower() and not i.D48_excess)]
        stds_TV03 = [i for i in analyses if (i.type == 'std' and 'tv03' in i.name.lower() and not i.D48_excess)]
        plt.figure(0)
        plt.figure(0).hold(True)
        plt.errorbar([CIT_Carrara_CRF for i in range(len(stds_Carrara))],[j.D47_CRF-CIT_Carrara_CRF for j in stds_Carrara], yerr = [k.D47_sterr for k in stds_Carrara], fmt = 'o')
        plt.errorbar([TV03_CRF for i in range(len(stds_TV03))],[j.D47_CRF-TV03_CRF for j in stds_TV03], yerr = [k.D47_sterr for k in stds_TV03], fmt = 'o')
        plt.xlim(0.3, 0.7)
        plt.xlabel(ur'$\mathrm{}\Delta_{47, CRF} \/ (\u2030)}$')
        plt.ylabel(ur'$\mathrm{\Delta_{47, measured}-\Delta_{47, expected} \/ (\u2030)}$')
        plt.show()

    return

def CI_hg_slope_finder(hgs):
    # hgs = [i for i in analyses if (i.type == 'hg' and not i.D48_excess)]
    d47_hgs = np.asarray([i.d47 for i in hgs])
    D47_raw_hgs = np.asarray([i.D47_raw for i in hgs])
    d47_stdev_hgs = np.asarray([i.d47_stdev for i in hgs])
    D47_sterr_hgs = np.asarray([i.D47_sterr for i in hgs])
    # Now, make hg D47 line, using a York regression
    global hg_slope, hg_intercept
    hg_slope, hg_intercept, r_47, sm, sb, xc, yc, ct = lsqcubic(d47_hgs, D47_raw_hgs,d47_stdev_hgs, D47_sterr_hgs)

    return hg_slope

#
# def CI_hg_slope_finder(d47_hgs, D47_raw_hgs, d47_stdev_hgs, D47_sterr_hgs):
#     # Now, make hg D47 line, using a York regression
#     global hg_slope, hg_intercept
#     hg_slope, hg_intercept, r_47, sm, sb, xc, yc, ct = lsqcubic(d47_hgs, D47_raw_hgs,d47_stdev_hgs, D47_sterr_hgs)
#
#     return hg_slope




def CI_48_excess_checker(analyses, showFigures = False):
    '''Checks for 48 excess using hg slope and int, based on a certain pre-specified tolerance'''
    D48_excess_tolerance = 1
    hgs = [i for i in analyses if i.type == 'hg']
    hg_slope_48, hg_intercept_48, r_48, sm_48, sb_48, xc_48, yc_, ct_ = lsqcubic(np.asarray([i.d48 for i in hgs]), np.asarray([i.D48_raw for i in hgs]),
    np.asarray([i.d48_stdev for i in hgs]), np.asarray([i.D48_stdev for i in hgs]))
    for i in analyses:
        D48_predicted = i.d48*hg_slope_48 + hg_intercept_48
        D48_excess_value = i.D48_raw - D48_predicted
        if np.abs(D48_excess_value) > 1:
            i.D48_excess = True
            print('D48 excess found for sample: ' + i.name + ', ' + str(i.num))


    if showFigures:
        plt.figure(1)
        plt.subplot(2,1,1)
        d48s = np.asarray([i.d48 for i in hgs])
        D48s_model = d48s * hg_slope_48 + hg_intercept_48
        # First, plotting the D48 line, TODO: include regression
        plt.figure(1).hold(True)
        plt.errorbar(np.asarray([i.d48 for i in hgs]), np.asarray([i.D48_raw for i in hgs]),
        xerr = np.asarray([i.d48_stdev for i in hgs]), yerr = np.asarray([i.D48_stdev for i in hgs]),
        fmt = 'bo')
        plt.plot(d48s, D48s_model,'r-')
        plt.xlabel(ur'$\delta^{48} \/ ( \u2030 $)')
        plt.ylabel(ur'$\Delta_{48} \/ ( \u2030 $)')

        # Plotting all samples in D47 vs d47 space,
        plt.subplot(2,1,2)
        plt.errorbar(np.asarray([i.d47 for i in analyses]), np.asarray([i.D47_raw for i in analyses]),
        xerr = np.asarray([i.d47_stdev for i in analyses]), yerr = np.asarray([i.D47_sterr for i in analyses]),
        fmt = 'o')
        plt.xlabel(ur'$\delta^{47} \/ ( \u2030 $)')
        plt.ylabel(ur'$\Delta_{47} \/ ( \u2030 $)')
        # plt.savefig('D48_line.pdf', format = 'pdf')
        plt.show()


    return




def CI_CRF_corrector(analysis, objName):
    '''Function to apply the heated gas correction in the Caltech Ref Frame'''

    acid_digestion_correction = 0.081
    try:
        D47_hg_corrected = analysis.D47_raw - (analysis.d47*hg_slope + hg_intercept)
        D47_stretching = D47_hg_corrected *(-0.8453)/hg_intercept
        D47_CRF = D47_stretching + acid_digestion_correction
    except NameError:
        return(np.NaN)

    return(D47_CRF)

def CI_hg_values(instance, objName):
    '''Placeholder for hg slope and int when used in CI classes'''
    try:
        values = {'hg_slope': hg_slope, 'hg_intercept': hg_intercept}
        return(values[objName])
    except(NameError, AttributeError):
        return(np.NaN)


def Daeron_data_creator(analyses):
    '''Creates a list of dictionaries in the format needed for M. Daeron's data
    processing script '''
    daeronData = []
    for i in analyses:
        daeronData.append({'label': i.name, 'd47': i.d47, 'rawD47': i.D47_raw,
        'srawD47': i.D47_sterr} )
        if i.type == 'hg':
            daeronData[-1]['TCO2eq'] = 1000.0
            daeronData[-1]['trueD47'] = xlcor47_modified.CO2eqD47(1000.0)
        elif i.type == 'eg':
            daeronData[-1]['TCO2eq'] = 25.0
            daeronData[-1]['trueD47'] = xlcor47_modified.CO2eqD47(25.0)
    return(daeronData)

def Daeron_data_processer(analyses, showFigures = False):
    '''Performs Daeron-style correction to put clumped isotope data into the ARF'''

    daeronData = Daeron_data_creator(analyses)

    global daeronBestFitParams, CorrelationMatrix
    daeronBestFitParams, CorrelationMatrix = xlcor47_modified.process_data(daeronData)

    for i in range(len(daeronData)):
        analyses[i].D47_ARF = daeronData[i]['corD47']
        analyses[i].D47_error_internal = daeronData[i]['scorD47_internal']
        analyses[i].D47_error_model = daeronData[i]['scorD47_model']
        analyses[i].D47_error_all = daeronData[i]['scorD47_all']

    if showFigures:
        xlcor47_modified.plot_data( daeronData, daeronBestFitParams, CorrelationMatrix ,'filename')

    return

def ExportSequence(analyses):
    '''Most common export sequence'''
    print('Exporting to temporary CIDS and FlatList files ')
    print('Exporting full acqs to a CIDS sheet...')
    exportNameCIDS = 'autoCIDS_Export'
    CIDS_exporter(analyses, exportNameCIDS)
    print('Exporting analyses to a flatlist...')
    exportNameFlatlist = 'autoFlatListExport'
    FlatList_exporter(analyses,exportNameFlatlist)
    print('Analyses successfully exported')
    doDaeron = raw_input('Export analyses for a Daeron-style ARF reduction (y/n)? ')
    if doDaeron.lower() == 'y':
        exportNameDaeron = 'autoDaeronExport'
        Get_gases(analyses)
        Daeron_exporter(analyses,exportNameDaeron)

    return

# def CI_background_correction(instance, objName):
#     voltSamTemp = np.copy(instance.voltSam_raw)
#     voltRefTemp = np.copy(instance.voltRef_raw)
#
#     slopeArray = np.array([0, 0 , 0, mass47PblSlope, 0, 0])
#     interceptArray  = np.array([0, 0, 0, mass47PblIntercept, 0, 0])
#
#     # mass 47 correction
#     try:
#         # print('Correcting voltages now')
#         # print('using this mass47 slope: '+ str(mass47PblSlope))
#         voltSamTemp[:,3] = voltSamTemp[:,3] - mass47PblSlope * voltSamTemp[:,0] - mass47PblIntercept
#         voltRefTemp[:,3] = voltRefTemp[:,3] - mass47PblSlope * voltRefTemp[:,0] - mass47PblIntercept
#     except(IndexError):
#         return 0
#
#     voltages = {'voltSam': voltSamTemp, 'voltRef': voltRefTemp}
#
#     return(voltages[objName])

def CI_background_correction(instance, objName):
    # voltSamTemp = np.copy(instance.voltSam_raw)
    # voltRefTemp = np.copy(instance.voltRef_raw)

    slopeArray = np.array([0, 0 , 0, mass47PblSlope, 0, 0])
    # interceptArray  = np.array([0, 0, 0, mass47PblIntercept, 0, 0])

    # mass 47 correction
    try:
        # print('Correcting voltages now')
        # print('using this mass47 slope: '+ str(mass47PblSlope))
        voltSamTemp = instance.voltSam_raw - np.transpose(np.tile(instance.voltSam_raw[:,0],(len(slopeArray),1)))*slopeArray
        voltRefTemp = instance.voltRef_raw - np.transpose(np.tile(instance.voltRef_raw[:,0],(len(slopeArray),1)))*slopeArray
    except(IndexError):
        return 0

    voltages = {'voltSam': voltSamTemp, 'voltRef': voltRefTemp}

    return(voltages[objName])


def Set_mass_47_pbl(pblSlope):
    global mass47PblSlope
    mass47PblSlope = pblSlope
    # mass47PblIntercept = pblIntercept




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


def lsqcubic(X, Y, sX, sY, tl=1e-6):
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
