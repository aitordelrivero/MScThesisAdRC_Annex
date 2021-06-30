#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

SCRIPT TO VISUALIZE THE CSV DATA

(Datasets available at https://wiki-bsse.ethz.ch/pages/viewpage.action?pageId=247308834)

"""

from SLiCAP import *
def Cadencecsv2traces(csvFile, absx = False, logx = False, absy = False, logy = False, limiter=None):
    """
    Generates a dictionary with traces (key = label, value = trace object) from
    data from a csv file. The CSV file should have the following structure:

    x0_label, y0_label, x1_label, y1_label, ...
    x0_0    , y0_0    , x1_0    , y1_0    , ...
    x0_1    , y0_1    , x1_1    , y1_1    , ...
    ...     , ...     , ...     , ...     , ...

    The traces will be named  with their y label.

    :param csvFile: name of the csv file (in the ini.csvPath directory)
    :type csvFile: str

    :param absx: if 'True', it applies the absolute (abs) function to the indpendent variable data (xData)
    :type bool

    :param logx: if 'True', it applies the logarithm in base 10 (log10) function to the independent variable data (xData)
    :type bool

    :param absy: if 'True', it applies the absolute (abs) function to the dependent variable data (yData)
    :type bool

    :param logy: if 'True', it applies the logarithm in base 10 (log10) function to the dependent variable data (yData)
    :type bool

    :return: dictionary with key-value pairs:

             - key: *str*: label of the trace
             - value: *SLiCAPplots.trace* trace object

    :rtype: dict
    """
    try:
        f = open(ini.csvPath + csvFile)
        lines = f.readlines()
        f.close()
    except:
        print('Error: could not find CSV trace data:', ini.csvPath + csvFile)
        return {}
    traceDict = {}
    labels = []
    for i in range(len(lines)):
        if i==0:
            data = lines[i][1:-1].split('","')
        else:
            data = lines[i].split(',')
        #print(data)
        if len(data) % 2 != 0:
            print("Error: expected an even number of columns in csv file:", ini.csvPath + csvFile)
            return traceDict
        elif i == 0:
            for j in range(int(len(data)/2)):
                labels.append(data[2*j+1])
                traceDict[data[2*j+1]] = trace([[], []])
                traceDict[data[2*j+1]].xData = []
                traceDict[data[2*j+1]].yData = []
        else:
            for j in range(len(labels)):
                try:
                    xData = eval(data[2*j])
                except:
                    xData = limiter #the limiter is the maximum value of the independent variable or None (by now it is given as a parameter)
                if absx:
                    xData = abs(xData)
                if logx:
                    try:
                        xData = np.log10(xData)
                    except:
                        print("Could not calculate the log10 of the xData of:", ini.csvPath + csvFile)
                try:
                    yData = eval(data[2*j+1])
                except:
                    yData = 0
                if absy:
                    yData = abs(yData)
                if logy:
                    try:
                        yData = np.log10(yData)
                    except:
                        print("Could not calculate the log10 of the yData of:", ini.csvPath + csvFile)
                traceDict[labels[j]].xData.append(xData)
                traceDict[labels[j]].yData.append(yData)
    #print("Available data")
    for label in labels:
            traceDict[label].xData = np.array(traceDict[label].xData)
            traceDict[label].yData = np.array(traceDict[label].yData)
            traceDict[label].label = label
    return traceDict

t1 = time()

prj = initProject('Visualization')


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#COMPARISON OF LOOP GAINS FOR DIFFERENT DESIGN STEPS. GAIN=200MV/A######################################################
select=True ###########################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("COMPARISON OF LOOP GAINS FOR DIFFERENT DESIGN STEPS GAIN 200MV A")

    Trace1 = Cadencecsv2traces('LoopGainControllerMag.csv')
    Trace2 = Cadencecsv2traces('LoopgainMagPrelayout3.csv')
    Trace3 = Cadencecsv2traces('LoopGainMagPrePostLayoutComparisonOpAmp.csv')
    Trace4 = Cadencecsv2traces('LoopGainMagPrePostLayoutComparisonWholeCircuit.csv')
    keys=list(Trace1.keys())
    #print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,Rfb=2e+08) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain dB20 (Vclamp=1.65,Rfb=2e+08) Y","Prelayout. Ideal resistor.")
            print(str(Trace1sel[key].label))
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("Loop Gain dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y","Prelayout. Pseudoresistor.")
            print(str(Trace2sel[(key+'x')].label))
    Trace3sel={}
    keys=list(Trace3.keys())
    for key in keys:
        if key in ['Loop Gain db20 Postlayout (Vclamp=1.65,Rfb=2e+08) Y','Loop Gain db20 Prelayout (Vclamp=1.65,Rfb=2e+08) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain db20 Postlayout (Vclamp=1.65,Rfb=2e+08) Y","Postlayout. Ideal resistor.")
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain db20 Prelayout (Vclamp=1.65,Rfb=2e+08) Y","Prelayout. Ideal resistor. Blackbox.")
            print(str(Trace3sel[key].label))
    Trace4sel={}
    keys=list(Trace4.keys())
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y','Loop Gain dB20 Postlayout (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain dB20 Postlayout (Vclamp=1.65,A=0,Idac=9e-08) Y","Postlayout. Whole circuit.")
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y","Prelayout. Pseudoresistor. Blackbox.")
            print(str(Trace4sel[key].label))
    Trace1sel.update(Trace2sel)
    Trace1sel.update(Trace3sel)
    Trace1sel.update(Trace4sel)
    print(list(Trace1sel.keys()))
    Graph = plot('loopgainscomparison_mag_200M', 'Loop gain magnitude (Gain=200MV/A)', 'semilogx', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Loop gain magnitude', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)


    Trace1 = Cadencecsv2traces('LoopGainControllerPhase.csv')
    Trace2 = Cadencecsv2traces('LoopgainPhasePrelayout3.csv')
    Trace3 = Cadencecsv2traces('LoopGainPhasePrePostLayoutComparisonOpAmp.csv')
    Trace4 = Cadencecsv2traces('LoopGainPhasePrePostLayoutComparisonWholeCircuit.csv')
    keys=list(Trace1.keys())
    #print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,Rfb=2e+08) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain Phase (Vclamp=1.65,Rfb=2e+08) Y","Prelayout. Ideal resistor.")
            print(str(Trace1sel[key].label))
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("Loop Gain Phase (Vclamp=1.65,A=0,Idac=9e-08) Y","Prelayout. Pseudoresistor.")
            print(str(Trace2sel[(key+'x')].label))
    Trace3sel={}
    keys=list(Trace3.keys())
    for key in keys:
        if key in ['Loop Gain Phase Postlayout (Vclamp=1.65,Rfb=2e+08) Y','Loop Gain Phase Prelayout (Vclamp=1.65,Rfb=2e+08) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain Phase Postlayout (Vclamp=1.65,Rfb=2e+08) Y","Postlayout. Ideal resistor.")
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain Phase Prelayout (Vclamp=1.65,Rfb=2e+08) Y","Prelayout. Ideal resistor. Blackbox.")
            print(str(Trace3sel[key].label))
    Trace4sel={}
    keys=list(Trace4.keys())
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,A=0,Idac=9e-08) Y','Loop Gain Phase Postlayout (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain Phase Postlayout (Vclamp=1.65,A=0,Idac=9e-08) Y","Postlayout. Whole circuit.")
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain Phase (Vclamp=1.65,A=0,Idac=9e-08) Y","Prelayout. Pseudoresistor. Blackbox.")
            print(str(Trace4sel[key].label))
    Trace1sel.update(Trace2sel)
    Trace1sel.update(Trace3sel)
    Trace1sel.update(Trace4sel)
    print(list(Trace1sel.keys()))
    Graph = plot('loopgainscomparison_phase_200M', 'Loop gain phase (Gain=200MV/A)', 'semilogx', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Loop gain phase', yUnits = '$deg$', xLim = [] , yLim = [], show = False)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#COMPARISON OF LOOP GAINS FOR DIFFERENT DESIGN STEPS. GAIN=2MV/A######################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("COMPARISON OF LOOP GAINS FOR DIFFERENT DESIGN STEPS GAIN 2 MV A")

    Trace1 = Cadencecsv2traces('LoopGainControllerMag.csv')
    Trace2 = Cadencecsv2traces('LoopgainMagPrelayout3.csv')
    Trace3 = Cadencecsv2traces('LoopGainMagPrePostLayoutComparisonOpAmp.csv')
    Trace4 = Cadencecsv2traces('LoopGainMagPrePostLayoutComparisonWholeCircuit.csv')
    keys=list(Trace1.keys())
    #print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,Rfb=2000000) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain dB20 (Vclamp=1.65,Rfb=2e+08) Y","Prelayout. Ideal resistor.")
            print(str(Trace1sel[key].label))
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y']:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("Loop Gain dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y","Prelayout. Pseudoresistor.")
            print(str(Trace2sel[(key+'x')].label))
    Trace3sel={}
    keys=list(Trace3.keys())
    for key in keys:
        if key in ['Loop Gain dB20 Postlayout (Vclamp=1.65,Rfb=2000000) Y','Loop Gain dB20 Prelayout (Vclamp=1.65,Rfb=2000000) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain dB20 Postlayout (Vclamp=1.65,Rfb=2000000) Y","Postlayout. Ideal resistor.")
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain dB20 Prelayout (Vclamp=1.65,Rfb=2000000) Y","Prelayout. Ideal resistor. Blackbox.")
            print(str(Trace3sel[key].label))
    Trace4sel={}
    keys=list(Trace4.keys())
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y','Loop Gain dB20 Postlayout (Vclamp=1.65,A=1,Idac=3.29e-07) Y']:
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y","Prelayout. Pseudoresistor. Blackbox.")
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain dB20 Postlayout (Vclamp=1.65,A=1,Idac=3.29e-07) Y","Postlayout. Whole circuit.")
            print(str(Trace4sel[key].label))
    Trace1sel.update(Trace2sel)
    Trace1sel.update(Trace3sel)
    Trace1sel.update(Trace4sel)
    print(list(Trace1sel.keys()))
    Graph = plot('loopgainscomparison_mag_2M', 'Loop gain magnitude (Gain=2MV/A)', 'semilogx', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Loop gain magnitude', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)


    Trace1 = Cadencecsv2traces('LoopGainControllerPhase.csv')
    Trace2 = Cadencecsv2traces('LoopgainPhasePrelayout3.csv')
    Trace3 = Cadencecsv2traces('LoopGainPhasePrePostLayoutComparisonOpAmp.csv')
    Trace4 = Cadencecsv2traces('LoopGainPhasePrePostLayoutComparisonWholeCircuit.csv')
    keys=list(Trace1.keys())
    #print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,Rfb=2000000) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain Phase (Vclamp=1.65,Rfb=2000000) Y","Prelayout. Ideal resistor.")
            print(str(Trace1sel[key].label))
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,A=1,Idac=3.29e-07) Y']:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("Loop Gain Phase (Vclamp=1.65,A=1,Idac=3.29e-07) Y","Prelayout. Pseudoresistor.")
            print(str(Trace2sel[(key+'x')].label))
    Trace3sel={}
    keys=list(Trace3.keys())
    for key in keys:
        if key in ['Loop Gain Phase Postlayout (Vclamp=1.65,Rfb=2000000) Y','Loop Gain Phase Prelayout (Vclamp=1.65,Rfb=2000000) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain Phase Postlayout (Vclamp=1.65,Rfb=2000000) Y","Postlayout. Ideal resistor.")
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain Phase Prelayout (Vclamp=1.65,Rfb=2000000) Y","Prelayout. Ideal resistor. Blackbox.")
            print(str(Trace3sel[key].label))
    Trace4sel={}
    keys=list(Trace4.keys())
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,A=1,Idac=3.29e-07) Y','Loop Gain Phase Postlayout (Vclamp=1.65,A=1,Idac=3.29e-07) Y']:
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain Phase (Vclamp=1.65,A=1,Idac=3.29e-07) Y","Prelayout. Pseudoresistor. Blackbox.")
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain Phase Postlayout (Vclamp=1.65,A=1,Idac=3.29e-07) Y","Postlayout. Whole circuit.")
            print(str(Trace4sel[key].label))
    Trace1sel.update(Trace2sel)
    Trace1sel.update(Trace3sel)
    Trace1sel.update(Trace4sel)
    print(list(Trace1sel.keys()))
    Graph = plot('loopgainscomparison_phase_2M', 'Loop gain phase (Gain=2MV/A)', 'semilogx', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Loop gain phase', yUnits = '$deg$', xLim = [] , yLim = [], show = False)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#COMPARISON OF LOOP GAINS FOR DIFFERENT DESIGN STEPS. GAIN=20G/A######################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("COMPARISON OF LOOP GAINS FOR DIFFERENT DESIGN STEPS GAIN 20 GV A")

    Trace1 = Cadencecsv2traces('LoopGainControllerMag.csv')
    Trace2 = Cadencecsv2traces('LoopgainMagPrelayout3.csv')
    Trace3 = Cadencecsv2traces('LoopGainMagPrePostLayoutComparisonOpAmp.csv')
    Trace4 = Cadencecsv2traces('LoopGainMagPrePostLayoutComparisonWholeCircuit.csv')
    keys=list(Trace1.keys())
    #print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,Rfb=2e+10) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain dB20 (Vclamp=1.65,Rfb=2e+10) Y","Prelayout. Ideal resistor.")
            print(str(Trace1sel[key].label))
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y']:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("Loop Gain dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y","Prelayout. Pseudoresistor.")
            print(str(Trace2sel[(key+'x')].label))
    Trace3sel={}
    keys=list(Trace3.keys())
    for key in keys:
        if key in ['Loop Gain db20 Postlayout (Vclamp=1.65,Rfb=2e+10) Y','Loop Gain db20 Prelayout (Vclamp=1.65,Rfb=2e+10) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain db20 Postlayout (Vclamp=1.65,Rfb=2e+10) Y","Postlayout. Ideal resistor.")
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain db20 Prelayout (Vclamp=1.65,Rfb=2e+10) Y","Prelayout. Ideal resistor. Blackbox.")
            print(str(Trace3sel[key].label))
    Trace4sel={}
    keys=list(Trace4.keys())
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y','Loop Gain dB20 Postlayout (Vclamp=1.65,A=0,Idac=1e-09) Y']:
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y","Prelayout. Pseudoresistor. Blackbox.")
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain dB20 Postlayout (Vclamp=1.65,A=0,Idac=1e-09) Y","Postlayout. Whole circuit.")
            print(str(Trace4sel[key].label))
    Trace1sel.update(Trace2sel)
    Trace1sel.update(Trace3sel)
    Trace1sel.update(Trace4sel)
    print(list(Trace1sel.keys()))
    Graph = plot('loopgainscomparison_mag_20G', 'Loop gain magnitude (Gain=20GV/A)', 'semilogx', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Loop gain magnitude', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)


    Trace1 = Cadencecsv2traces('LoopGainControllerPhase.csv')
    Trace2 = Cadencecsv2traces('LoopgainPhasePrelayout3.csv')
    Trace3 = Cadencecsv2traces('LoopGainPhasePrePostLayoutComparisonOpAmp.csv')
    Trace4 = Cadencecsv2traces('LoopGainPhasePrePostLayoutComparisonWholeCircuit.csv')
    keys=list(Trace1.keys())
    #print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,Rfb=2e+10) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain Phase (Vclamp=1.65,Rfb=2e+10) Y","Prelayout. Ideal resistor.")
            print(str(Trace1sel[key].label))
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,A=0,Idac=1e-09) Y']:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("Loop Gain Phase (Vclamp=1.65,A=0,Idac=1e-09) Y","Prelayout. Pseudoresistor.")
            print(str(Trace2sel[(key+'x')].label))
    Trace3sel={}
    keys=list(Trace3.keys())
    for key in keys:
        if key in ['Loop Gain Phase Postlayout (Vclamp=1.65,Rfb=2e+10) Y','Loop Gain Phase Prelayout (Vclamp=1.65,Rfb=2e+10) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain Phase Postlayout (Vclamp=1.65,Rfb=2e+10) Y","Postlayout. Ideal resistor.")
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain Phase Prelayout (Vclamp=1.65,Rfb=2e+10) Y","Prelayout. Ideal resistor. Blackbox.")
            print(str(Trace3sel[key].label))
    Trace4sel={}
    keys=list(Trace4.keys())
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,A=0,Idac=1e-09) Y','Loop Gain Phase Postlayout (Vclamp=1.65,A=0,Idac=1e-09) Y']:
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain Phase (Vclamp=1.65,A=0,Idac=1e-09) Y","Prelayout. Pseudoresistor. Blackbox.")
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain Phase Postlayout (Vclamp=1.65,A=0,Idac=1e-09) Y","Postlayout. Whole circuit.")
            print(str(Trace4sel[key].label))
    Trace1sel.update(Trace2sel)
    Trace1sel.update(Trace3sel)
    Trace1sel.update(Trace4sel)
    print(list(Trace1sel.keys()))
    Graph = plot('loopgainscomparison_phase_20G', 'Loop gain phase (Gain=20GV/A)', 'semilogx', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Loop gain phase', yUnits = '$deg$', xLim = [] , yLim = [], show = False)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)


    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    #NOISE ######################################################
select=True#############################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
if select:
    htmlPage("Noise comparison")
    Trace1 = Cadencecsv2traces('NoiseController.csv')
    Trace2 = Cadencecsv2traces('NoiseSpectrumPrelayout3.csv')
    Trace3 = Cadencecsv2traces('NoiseSpectrumPostlayout.csv')
    keys=list(Trace1.keys())
    #print(str(keys))
    Trace1sel={}
    # for key in keys:
    #     if key in ['onoise (Vclamp=1.65,Rfb=2e+08) Y']:
    #         Trace1sel[key]=Trace1[key]
    #         Trace1sel[key].label=Trace1sel[key].label.replace("onoise (Vclamp=1.65,Rfb=2000000) Y","Controller - Gain 2M")
    #         Trace1sel[key].label=Trace1sel[key].label.replace("onoise (Vclamp=1.65,Rfb=2e+08) Y","Controller - Gain 200M")
    #         Trace1sel[key].label=Trace1sel[key].label.replace("onoise (Vclamp=1.65,Rfb=2e+10) Y","Controller - Gain 20G")
    #         print(str(Trace1sel[key].label))
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['onoise (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("onoise (Vclamp=1.65,A=1,Idac=3.27e-7) Y","Total noise, Prelayout. - Gain 2M")
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("onoise (Vclamp=1.65,A=0,Idac=9e-08) Y","Total noise, Prelayout. - Gain 200M")
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("onoise (Vclamp=1.65,A=0,Idac=1e-9) Y","Total noise, Prelayout. - Gain 20G")
            print(str(Trace2sel[(key+'x')].label))
    Trace1sel.update(Trace2sel)
    Trace3sel={}
    keys=list(Trace3.keys())
    for key in keys:
        if key in ['onoise (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace3sel[(key+'y')]=Trace3[key]
            Trace3sel[(key+'y')].label=Trace3sel[(key+'y')].label.replace("onoise (Vclamp=1.65,A=1,Idac=3.29e-07) Y","Total noise, Postlayout. - Gain 2M")
            Trace3sel[(key+'y')].label=Trace3sel[(key+'y')].label.replace("onoise (Vclamp=1.65,A=0,Idac=9e-08) Y","Total noise, Postlayout. - Gain 200M")
            Trace3sel[(key+'y')].label=Trace3sel[(key+'y')].label.replace("onoise (Vclamp=1.65,A=0,Idac=1e-09) Y","Total noise, Postlayout. - Gain 20G")
            print(str(Trace3sel[(key+'y')].label))
    Trace1sel.update(Trace3sel)
    print(list(Trace1sel.keys()))
    Graph = plot('Noisecomparison', 'Onoise', 'log', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Noise', yUnits = '$V**2/Hz$', xLim = [] , yLim = [], show = False)
    #del Trace1, Trace2, Trace3, Trace1sel, Trace2sel, Trace3sel
    fig2html(Graph, 600)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#NOISE RMS #############################################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("NOISE RMS COMPARISON")
    Trace1 = csv2traces('NoiseRms3.csv')
    Trace2 = csv2traces('RMSNoiseController.csv')
    Trace3 = csv2traces('RMSNoisePostLayout.csv')
    keys=list(Trace1.keys())
    print(keys)
    Trace1sel={}
    for key in keys:
        if key in ['leafValue( Vrms_onoise "A" 1 ) vs leafValue( ReqAC_Prelayout "A" 1 ) Y','leafValue( Vrms_onoise "A" 0 ) vs leafValue( ReqAC_Prelayout "A" 0 ) Y\n']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace('leafValue( Vrms_onoise "A" 0 ) vs leafValue( ReqAC_Prelayout "A" 0 ) Y\n',"TIA High gain mode - Prelayout")
            Trace1sel[key].label=Trace1sel[key].label.replace('leafValue( Vrms_onoise "A" 1 ) vs leafValue( ReqAC_Prelayout "A" 1 ) Y',"TIA Low gain mode - Prelayout")
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ['Vrms_onoise (Vclamp=1.65) Y']:
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("Vrms_onoise (Vclamp=1.65) Y","Controller - Prelayout")
            print(str(Trace2sel[key].label))
    Trace1sel.update(Trace2sel)
    keys=list(Trace3.keys())
    Trace3sel={}
    for key in keys:
        if key in ['leafValue( Vrms_onoise "A" 1 ) vs leafValue( Req "A" 1 ) Y\n','leafValue( Vrms_onoise "A" 0 ) vs leafValue( Req "A" 0 ) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace('leafValue( Vrms_onoise "A" 1 ) vs leafValue( Req "A" 1 ) Y',"TIA Low gain mode - Postlayout")
            Trace3sel[key].label=Trace3sel[key].label.replace('leafValue( Vrms_onoise "A" 0 ) vs leafValue( Req "A" 0 ) Y',"TIA High gain mode - Postlayout")
            print(str(Trace3sel[key].label))
    Trace1sel.update(Trace3sel)
    GraphNoiseRMS = plot('NoiseRMScomparison', 'Output referred noise - RMS', 'log', Trace1sel, xName = 'Gain', xUnits = 'V/A', yName = 'onoise', yUnits = 'Vrms',xLim = [] , yLim = [], show = False)
    fig2html(GraphNoiseRMS, 600)
    #del Trace1, Trace2, Trace3, Trace1sel, Trace2sel, Trace3sel


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#TRANSIENTS 200M#############################################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("COMBINED TRANSIENTS GAIN 200MV A")
    Trace1 = Cadencecsv2traces('MidgainTransientPrelayout3.csv',limiter=0.2)
    Trace2 = Cadencecsv2traces('ControllerTransient.csv',limiter=0.2)
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['/out (Vclamp_delta=0.375,tr=0.1) Y','/out (Vclamp_delta=0,tr=0.01) Y','/out (Vclamp_delta=1.5,tr=0.1) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=0.375,tr=0.1) Y","Pseudoresistor. Prelayout. Vclamp_delta=0.375,tr=0.1")
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=0,tr=0.01) Y","Pseudoresistor. Prelayout. Vclamp_delta=0,tr=0.01")
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=1.5,tr=0.1) Y","Pseudoresistor. Prelayout. Vclamp_delta=1.5,tr=0.1")
    keys=list(Trace2.keys())
    #print(str(keys))
    Trace2sel={}
    for key in keys:
        if key in ['/out (Im_pk=5e-09,Vclamp_delta=0.375,tr=0.1) Y','/out (Im_pk=5e-09,Vclamp_delta=0,tr=0.01) Y','/out (Im_pk=5e-09,Vclamp_delta=1.5,tr=0.1) Y']:
            #print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("/out (Im_pk=5e-09,Vclamp_delta=0.375,tr=0.1) Y","Ideal resistor. Prelayout. Vclamp_delta=0.375,tr=0.1")
            Trace2sel[key].label=Trace2sel[key].label.replace("/out (Im_pk=5e-09,Vclamp_delta=0,tr=0.01) Y","Ideal resistor. Prelayout. Vclamp_delta=0,tr=0.01")
            Trace2sel[key].label=Trace2sel[key].label.replace("/out (Im_pk=5e-09,Vclamp_delta=1.5,tr=0.1) Y","Ideal resistor. Prelayout. Vclamp_delta=1.5,tr=0.1")
    Trace1sel.update(Trace2sel)
    Graph = plot('transients200M', 'Output signal amplitude (Gain=200MV/A)', 'lin', Trace1sel, xName = 't', xUnits = 's', yName = 'Amplitude', yUnits = '$V$', xLim = [0, 0.199] , yLim = [], show = False)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#TRANSIENTS 2M########################################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("COMBINED TRANSIENTS GAIN 2MV A")
    Trace1 = Cadencecsv2traces('LowgainTransientPrelayout3.csv',limiter=0.2)
    Trace2 = Cadencecsv2traces('ControllerTransient.csv',limiter=0.2)
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['/out (Vclamp_delta=0.375,tr=0.001) Y','/out (Vclamp_delta=0.75,tr=0.01) Y','/out (Vclamp_delta=0,tr=0.1) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=0.375,tr=0.001) Y","Pseudoresistor. Prelayout. Vclamp_delta=0.375,tr=0.001")
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=0.75,tr=0.01) Y","Pseudoresistor. Prelayout. Vclamp_delta=0.75,tr=0.01")
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=0,tr=0.1) Y","Pseudoresistor. Prelayout. Vclamp_delta=0,tr=0.1")
    keys=list(Trace2.keys())
    #print(str(keys))
    Trace2sel={}
    for key in keys:
        if key in ['/out (Im_pk=5e-07,Vclamp_delta=0.375,tr=0.001) Y','/out (Im_pk=5e-07,Vclamp_delta=0.75,tr=0.01) Y','/out (Im_pk=5e-07,Vclamp_delta=0,tr=0.1) Y']:
            #print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("/out (Im_pk=5e-07,Vclamp_delta=0.375,tr=0.001) Y","Ideal resistor. Prelayout. Vclamp_delta=0.375,tr=0.001")
            Trace2sel[key].label=Trace2sel[key].label.replace("/out (Im_pk=5e-07,Vclamp_delta=0.75,tr=0.01) Y","Ideal resistor. Vclamp_delta=0.75,tr=0.01")
            Trace2sel[key].label=Trace2sel[key].label.replace("/out (Im_pk=5e-07,Vclamp_delta=0,tr=0.1) Y","Ideal resistor. Vclamp_delta=0,tr=0.1")
    Trace1sel.update(Trace2sel)
    Graph = plot('transients2M', 'Output signal amplitude (Gain=2MV/A)', 'lin', Trace1sel, xName = 't', xUnits = 's', yName = 'Amplitude', yUnits = '$V$', xLim = [0, 0.199] , yLim = [], show = False)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#TRANSIENTS 20G#########################################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("COMBINED TRANSIENTS GAIN 20GV A")
    Trace1 = Cadencecsv2traces('HighgainTransientPrelayout3.csv')
    Trace2 = Cadencecsv2traces('ControllerTransient.csv')
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['/out (Vclamp_delta=0.375,tr=0.001) Y','/out (Vclamp_delta=1.125,tr=1) Y','/out (Vclamp_delta=0,tr=0.1) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=0.375,tr=0.001) Y","Pseudoresistor. Prelayout. Vclamp_delta=0.375,tr=0.001")
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=1.125,tr=1) Y","Pseudoresistor. Prelayout. Vclamp_delta=1.125,tr=1")
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=0,tr=0.1) Y","Pseudoresistor. Prelayout. Vclamp_delta=0,tr=0.1")
    keys=list(Trace2.keys())
    #print(str(keys))
    Trace2sel={}
    for key in keys:
        if key in ['/out (Im_pk=5e-11,Vclamp_delta=0.375,tr=0.001) Y','/out (Im_pk=5e-11,Vclamp_delta=1.125,tr=1) Y','/out (Im_pk=5e-11,Vclamp_delta=0,tr=0.1) Y']:
            #print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("/out (Im_pk=5e-07,Vclamp_delta=0.375,tr=0.001) Y","Ideal resistor. Prelayout. Vclamp_delta=0.375,tr=0.001")
            Trace2sel[key].label=Trace2sel[key].label.replace("/out (Im_pk=5e-07,Vclamp_delta=0.75,tr=0.01) Y","Ideal resistor. Vclamp_delta=0.75,tr=0.01")
            Trace2sel[key].label=Trace2sel[key].label.replace("/out (Im_pk=5e-07,Vclamp_delta=0,tr=0.1) Y","Ideal resistor. Vclamp_delta=0,tr=0.1")
    Trace1sel.update(Trace2sel)
    Graph = plot('transients20G', 'Output signal amplitude (Gain=20GV/A)', 'lin', Trace1sel, xName = 't', xUnits = 's', yName = 'Amplitude', yUnits = '$V$', xLim = [0, 0.199] , yLim = [], show = False)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#INSTANT VOLTAGE CLAMP #########################################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("INSTANT VOLTAGE CLAMP")
    Trace1 = Cadencecsv2traces('TransientVoltageClampReferenceClamp_delta02_st_PreLayout.csv', limiter=0.0001)
    Trace2 = Cadencecsv2traces('TransientVoltageClampShortTao.csv', limiter=0.0001)
    keys=list(Trace1.keys())
    #print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['/inA (A=0,Idac=9.501601e-08,Vclamp=1.65) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("/inA (A=0,Idac=9.501601e-08,Vclamp=1.65) Y","Fast voltage clamp - TIA prelayout")
    keys=list(Trace2.keys())
    #print(str(keys))
    Trace2sel={}
    for key in keys:
        if key in ['/net8 (Rfb=2e+08,Vclamp=1.65) Y']:
            #print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("/net8 (Rfb=2e+08,Vclamp=1.65) Y","Fast voltage clamp - Controller prelayout")
    Trace1sel.update(Trace2sel)
    Graph = plot('fastvoltageclamp', 'Fast clamp (Gain=200MV/A)', 'lin', Trace1sel, xName = 't', xUnits = 's', yName = 'Amplitude', yUnits = '$V$', xLim = [0,0.0000999] , yLim = [], show = False)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#AC #########################################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("AC TRANSFER")
    Trace1 = Cadencecsv2traces('GainControllerMag.csv') #Ideal resistor
    Trace2 = Cadencecsv2traces('GainControllerMagPrelayout.csv')
    Trace3 = Cadencecsv2traces('GainControllerMagPostlayout.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['Gain dB20 (Vclamp=1.65,Rfb=2e+08) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Gain dB20 (Vclamp=1.65,Rfb=2e+08) Y","AC Mag Vclamp=1.65 Rfb=200M - Controller prelayout")
    Graph = plot('actransfer', 'AC transfer (Gain=200MV/A)', 'semilogx', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'AC', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ['']:
            print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("","")
    Trace1sel.update(Trace2sel)
    Graph = plot('actransfer', 'AC transfer (Gain=200MV/A)', 'semilogx', Trace2sel, xName = 'f', xUnits = 'Hz', yName = 'AC', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel



########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#PSRR summary table  #########################################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("PSRRi")
    #Trace1 = Cadencecsv2traces('PSRRv.csv')
    Trace4 = Cadencecsv2traces('PSRRiPreLayout.csv')
    Trace3 = Cadencecsv2traces('PSRRiPostLayout.csv')
    #keys=list(Trace1.keys())
    #print(str(keys))
    # Trace1sel={}
    # for key in keys:
    #     if key in ['PSRRv (Vclamp=1.65,Rfb=2e+08) Y']:
    #         print(key)
    #         Trace1sel[key]=Trace1[key]
    #         Trace1sel[key].label=Trace1sel[key].label.replace("PSRRv (Vclamp=1.65,Rfb=2e+08) Y","PSRRv Vclamp=1.65 Rfb=200M - Controller prelayout")
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['PSRRi (Vclamp=1.65,A=0,Idac=9e-08) Y','PSRRi (Vclamp=1.65,A=0,Idac=3.29e-07) Y','PSRRi (Vclamp=1.65,A=0,Idac=1e-09) Y','PSRRi (Vclamp=1.65,A=1,Idac=9e-08) Y','PSRRi (Vclamp=1.65,A=1,Idac=3.29e-07) Y','PSRRi (Vclamp=1.65,A=1,Idac=1e-09) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("PSRRi (Vclamp=1.65,A=0,Idac=9e-08) Y","Idac=9e-08, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("PSRRi (Vclamp=1.65,A=0,Idac=3.29e-07) Y","Idac=3.27e-07, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('PSRRi (Vclamp=1.65,A=0,Idac=1e-09) Y',"Idac=1e-09, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('PSRRi (Vclamp=1.65,A=1,Idac=9e-08) Y',"Idac=9e-08, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('PSRRi (Vclamp=1.65,A=1,Idac=3.29e-07) Y',"Idac=3.27e-07, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('PSRRi (Vclamp=1.65,A=1,Idac=1e-09) Y',"Idac=1e-09, LR mode")
    Graph = plot('PSRRi_prelayout', 'PSRRi - Pre layout', 'semilogx', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'PSRRv', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    keys=list(Trace3.keys())
    Trace3sel={}
    for key in keys:
        if key in ['PSRRi (Vclamp=1.65,A=0,Idac=9e-08) Y','PSRRi (Vclamp=1.65,A=0,Idac=3.29e-07) Y','PSRRi (Vclamp=1.65,A=0,Idac=1e-09) Y','PSRRi (Vclamp=1.65,A=1,Idac=9e-08) Y','PSRRi (Vclamp=1.65,A=1,Idac=3.29e-07) Y','PSRRi (Vclamp=1.65,A=1,Idac=1e-09) Y']:
            print(key)
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("PSRRi (Vclamp=1.65,A=0,Idac=9e-08) Y","Idac=9e-08, HR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace("PSRRi (Vclamp=1.65,A=0,Idac=3.29e-07) Y","Idac=3.27e-07, HR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace('PSRRi (Vclamp=1.65,A=0,Idac=1e-09) Y',"Idac=1e-09, HR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace('PSRRi (Vclamp=1.65,A=1,Idac=9e-08) Y',"Idac=9e-08, LR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace('PSRRi (Vclamp=1.65,A=1,Idac=3.29e-07) Y',"Idac=3.27e-07, LR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace('PSRRi (Vclamp=1.65,A=1,Idac=1e-09) Y',"Idac=1e-09, LR mode")
    Graph = plot('PSRRi_postlayout', 'PSRRi - Post layout', 'semilogx', Trace3sel, xName = 'f', xUnits = 'Hz', yName = 'PSRRv', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)

if select:
    htmlPage("PSRRv")
    #Trace1 = Cadencecsv2traces('PSRRv.csv')
    Trace4 = Cadencecsv2traces('PSRRvPreLayout.csv')
    Trace3 = Cadencecsv2traces('PSRRvPostLayout.csv')
    #keys=list(Trace1.keys())
    #print(str(keys))
    # Trace1sel={}
    # for key in keys:
    #     if key in ['PSRRv (Vclamp=1.65,Rfb=2e+08) Y']:
    #         print(key)
    #         Trace1sel[key]=Trace1[key]
    #         Trace1sel[key].label=Trace1sel[key].label.replace("PSRRv (Vclamp=1.65,Rfb=2e+08) Y","PSRRv Vclamp=1.65 Rfb=200M - Controller prelayout")
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['PSRRv (Vclamp=1.65,A=0,Idac=9e-08) Y','PSRRv (Vclamp=1.65,A=0,Idac=3.29e-07) Y','PSRRv (Vclamp=1.65,A=0,Idac=1e-09) Y','PSRRv (Vclamp=1.65,A=1,Idac=9e-08) Y','PSRRv (Vclamp=1.65,A=1,Idac=3.29e-07) Y','PSRRv (Vclamp=1.65,A=1,Idac=1e-09) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("PSRRv (Vclamp=1.65,A=0,Idac=9e-08) Y","Idac=9e-08, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("PSRRv (Vclamp=1.65,A=0,Idac=3.29e-07) Y","Idac=3.27e-07, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('PSRRv (Vclamp=1.65,A=0,Idac=1e-09) Y',"Idac=1e-09, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('PSRRv (Vclamp=1.65,A=1,Idac=9e-08) Y',"Idac=9e-08, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('PSRRv (Vclamp=1.65,A=1,Idac=3.29e-07) Y',"Idac=3.27e-07, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('PSRRv (Vclamp=1.65,A=1,Idac=1e-09) Y',"Idac=1e-09, LR mode")
    Graph = plot('PSRRv_prelayout', 'PSRRv - Pre layout', 'semilogx', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'PSRRv', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    keys=list(Trace3.keys())
    Trace3sel={}
    for key in keys:
        if key in ['PSRRv (Vclamp=1.65,A=0,Idac=9e-08) Y','PSRRv (Vclamp=1.65,A=0,Idac=3.29e-07) Y','PSRRv (Vclamp=1.65,A=0,Idac=1e-09) Y','PSRRv (Vclamp=1.65,A=1,Idac=9e-08) Y','PSRRv (Vclamp=1.65,A=1,Idac=3.29e-07) Y','PSRRv (Vclamp=1.65,A=1,Idac=1e-09) Y']:
            print(key)
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("PSRRv (Vclamp=1.65,A=0,Idac=9e-08) Y","Idac=9e-08, HR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace("PSRRv (Vclamp=1.65,A=0,Idac=3.29e-07) Y","Idac=3.27e-07, HR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace('PSRRv (Vclamp=1.65,A=0,Idac=1e-09) Y',"Idac=1e-09, HR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace('PSRRv (Vclamp=1.65,A=1,Idac=9e-08) Y',"Idac=9e-08, LR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace('PSRRv (Vclamp=1.65,A=1,Idac=3.29e-07) Y',"Idac=3.27e-07, LR mode")
            Trace3sel[key].label=Trace3sel[key].label.replace('PSRRv (Vclamp=1.65,A=1,Idac=1e-09) Y',"Idac=1e-09, LR mode")
    Graph = plot('PSRRv_postlayout', 'PSRRv - Post layout', 'semilogx', Trace3sel, xName = 'f', xUnits = 'Hz', yName = 'PSRRv', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#PSRRi #########################################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("PSRRi")
    Trace1 = Cadencecsv2traces('PSRRi.csv')
    Trace2 = Cadencecsv2traces('PSRRiPrelayout.csv')
    Trace3 = Cadencecsv2traces('PSRRiPostlayout.csv')
    keys=list(Trace1.keys())
    #print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['PSRRi (Vclamp=1.65,Rfb=2e+08) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("PSRRi (Vclamp=1.65,Rfb=2e+08) Y","PSRRi Vclamp=1.65 Rfb=200M - Controller prelayout")
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ['']:
            print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("","")
    Trace1sel.update(Trace2sel)
    Graph = plot('PSRRi', 'PSRRi (Gain=200MV/A)', 'semilogx', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'PSRRi', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel
    fig2html(Graph, 600)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel



########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#DC PSEUDORESISTANCE Theoretical Vc ####################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if select:
    htmlPage("DC pseudoresistance 8 vs 4 blocks")
    Trace1 = csv2traces('pseudoresistanceIdeal8Blocks.csv')
    Trace2 = csv2traces('pseudoresistanceIdeal4Blocks.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['ExperimentCompounREq (Vc=0.897) Y','ExperimentCompounREq (Vc=0.273) Y']:
            #print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("ExperimentCompounREq (", " ").replace(") Y","")+' 8 blocks'
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['ExperimentCompounREq (Vc=0.246) Y','ExperimentCompounREq (Vc=0.75) Y']:
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("ExperimentCompounREq (", " ").replace(") Y","")+' 4 blocks'
            #print(key)
    Trace1sel.update(Trace2sel) #merge dictionaries
    print(list(Trace1sel.keys()))
    Graph = plot('pseudoresistanceIdeal', 'Comparison between 4 blocks and 8 blocks', 'semilogy', Trace1sel, xName = 'Vin', xUnits = 'V', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [-2,2] , yLim = [], show = False)
    fig2html(Graph, 600)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#AC PSEUDORESISTANCE Theoretical Vc#####################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if select:
    htmlPage("AC pseudoresistance")
    Trace1 = csv2traces('pseudoresistanceIdeal8BlockAC.csv')
    Trace2 = csv2traces('pseudoresistanceIdeal4BlockAC.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['ReqAC (Vc=0.897) Y','ReqAC (Vc=0.273) Y']:
            #print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("ReqAC (", " ").replace(") Y","")+' 8 blocks'
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['ReqAC (Vc=0.246) Y','ReqAC (Vc=0.75) Y']:
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("ReqAC (", " ").replace(") Y","")+' 4 blocks'
            #print(key)
    Trace1sel.update(Trace2sel) #merge dictionaries
    print(list(Trace1sel.keys()))
    Graph = plot('pseudoresistanceIdealAC', 'Comparison between 4 blocks and 8 blocks', 'log', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#DC PSEUDORESISTANCE comparison     ####################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if select:
    htmlPage("DC pseudoresistance Gain 245M")
    Trace1 = csv2traces('pseudoresistanceIdeal8Blocks.csv') #Ideal Vc
    Trace2 = Cadencecsv2traces('pseudoresistanceIdealBias.csv')   #Ideal Ibias
    Trace3 = Cadencecsv2traces('ReqFinalDC_full.csv') #With DAC
    Trace4 = Cadencecsv2traces('ReqStandardDC_extension.csv') #Replacing the DAC by an ideal source (standard measurements)
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['ExperimentCompounREq (Vc=0.45) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("ExperimentCompounREq (Vc=0.45) Y","Ideal Vc")
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ['Req_BothBias (Vclamp=1.65,A=0,Ic=7.993965e-08) Y']:
            print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("Req_BothBias (Vclamp=1.65,A=0,Ic=7.993965e-08) Y","Ideal Ibias")
    Trace1sel.update(Trace2sel)
    keys=list(Trace3.keys())
    Trace3sel={}
    for key in keys:
        if key in ['Req (Vclamp=1.65,A=0,B=0,C=0,D=1,E=1,F=0,G=1,H=1) Y']:
            print(key)
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("Req (Vclamp=1.65,A=0,B=0,C=0,D=1,E=1,F=0,G=1,H=1) Y","With DAC")
    Trace1sel.update(Trace3sel)
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['Req (Vclamp=1.65,A=0,Ic=7.296614e-08) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Req (Vclamp=1.65,A=0,Ic=7.296614e-08) Y","Replacing DAC with current source (standard)")
    Trace1sel.update(Trace4sel)
    Graph = plot('pseudoresistanceIdealDC', 'Comparison DC response pseudoresistor', 'semilogy', Trace1sel, xName = 'Vin', xUnits = 'V', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#AC PSEUDORESISTANCE comparison    #####################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if select:
    htmlPage("AC pseudoresistance Gain 245M")
    Trace1 = csv2traces('pseudoresistanceIdeal8BlockAC.csv') #Ideal Vc
    Trace2 = Cadencecsv2traces('pseudoresistanceIdealBiasAC.csv')   #Ideal Ibias
    Trace3 = Cadencecsv2traces('ReqFinalAC_full.csv') #With DAC
    Trace4 = Cadencecsv2traces('ReqStandardAC_extension.csv') #Replacing the DAC by an ideal source (standard measurements)
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['ReqAC (Vc=0.45) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("ReqAC (Vc=0.45) Y","Ideal Vc")
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ['ReqAC (Vclamp=1.65,A=0,Ic=7.993965e-08) Y']:
            print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("ReqAC (Vclamp=1.65,A=0,Ic=7.993965e-08) Y","Ideal Ibias")
    Trace1sel.update(Trace2sel)
    keys=list(Trace3.keys())
    Trace3sel={}
    for key in keys:
        if key in ['ReqAC (Vclamp=1.65,A=0,B=0,C=0,D=1,E=1,F=0,G=1,H=1) Y']:
            print(key)
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("ReqAC (Vclamp=1.65,A=0,B=0,C=0,D=1,E=1,F=0,G=1,H=1) Y","With DAC")
    Trace1sel.update(Trace3sel)
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['ReqAC (Vclamp=1.65,A=0,Ic=7.296614e-08) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("ReqAC (Vclamp=1.65,A=0,Ic=7.296614e-08) Y","Replacing DAC with current source (standard)")
    Trace1sel.update(Trace4sel)
    Graph = plot('pseudoresistanceIdealAC', 'Comparison AC response pseudoresistor', 'log', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#AC PSEUDORESISTANCE summary table  #####################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if select:
    htmlPage("AC pseudoresistance")
    Trace4 = Cadencecsv2traces('ReqStandardAC_extension.csv') #Replacing the DAC by an ideal source (standard measurements)
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['ReqAC (Vclamp=1.65,A=0,Ic=7.296614e-08) Y', 'ReqAC (Vclamp=1.65,A=1,Ic=3.29e-07) Y', 'ReqAC (Vclamp=1.65,A=1,Ic=7.296614e-08) Y','ReqAC (Vclamp=1.65,A=1,Ic=1e-09) Y','ReqAC (Vclamp=1.65,A=0,Ic=3.29e-07) Y','ReqAC (Vclamp=1.65,A=0,Ic=1e-09) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("ReqAC (Vclamp=1.65,A=0,Ic=7.296614e-08) Y","Gain 245MV/A, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("ReqAC (Vclamp=1.65,A=1,Ic=3.29e-07) Y","Gain 1.85MV/A, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('ReqAC (Vclamp=1.65,A=1,Ic=7.296614e-08) Y',"Gain 4.03MV/A, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('ReqAC (Vclamp=1.65,A=1,Ic=1e-09) Y',"Gain 165MV/A, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('ReqAC (Vclamp=1.65,A=0,Ic=3.29e-07) Y',"Gain 64MV/A, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('ReqAC (Vclamp=1.65,A=0,Ic=1e-09) Y',"Gain 18.3GV/A, HR mode")
    Graph = plot('pseudoresistanceAC_summarytable', 'Pseudo-resistance AC response, Vclamp=1.65', 'log', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    del Trace4, Trace4sel

if select:
    Trace4 = Cadencecsv2traces('PostLayoutPseudoResistorAC.csv') #Replacing the DAC by an ideal source (standard measurements)
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['ReqAC (Vclamp=1.65,Ic=7.296614e-08,A=0) Y', 'ReqAC (Vclamp=1.65,Ic=3.29e-07,A=1) Y', 'ReqAC (Vclamp=1.65,Ic=7.296614e-08,A=1) Y','ReqAC (Vclamp=1.65,Ic=1e-09,A=1) Y','ReqAC (Vclamp=1.65,Ic=3.29e-07,A=0) Y','ReqAC (Vclamp=1,Ic=1e-09,A=0) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("ReqAC (Vclamp=1.65,Ic=7.296614e-08,A=0) Y","Gain 245MV/A, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("ReqAC (Vclamp=1.65,Ic=3.29e-07,A=1) Y","Gain 1.85MV/A, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('ReqAC (Vclamp=1.65,Ic=7.296614e-08,A=1) Y',"Gain 4.03MV/A, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('ReqAC (Vclamp=1.65,Ic=1e-09,A=1) Y',"Gain 165MV/A, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('ReqAC (Vclamp=1.65,Ic=3.29e-07,A=0) Y',"Gain 64MV/A, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('ReqAC (Vclamp=1,Ic=1e-09,A=0) Y',"Gain 18.3GV/A, HR mode")
    Graph = plot('pseudoresistanceAC_summarytable2', 'Pseudo-resistance AC response, Vclamp=1.65', 'log', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    del Trace4, Trace4sel


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#SERVOS summary table ##################################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if select:
    htmlPage("Servo s")
    Trace4 = Cadencecsv2traces('ServoMagPrePostLayoutComparisonWholeCircuit.csv') #Replacing the DAC by an ideal source (standard measurements)
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['Servo dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=3.29e-07) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=9e-08) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=1e-09) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y","Idac=9e-08, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=3.29e-07) Y","Idac=3.27e-07, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y',"Idac=1e-09, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=9e-08) Y',"Idac=9e-08, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y',"Idac=3.27e-07, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=1e-09) Y',"Idac=1e-09, LR mode")
    Graph = plot('servo_postlayout', 'Servo function - Post layout, Vclamp=1.65', 'semilogx', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'Servo Magnitude', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    del Trace4, Trace4sel
    Trace4 = Cadencecsv2traces('ServoMagPrelayout3.csv') #Replacing the DAC by an ideal source (standard measurements)
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['Servo dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=3.29e-07) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=9e-08) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=1e-09) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y","Idac=9e-08, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=3.29e-07) Y","Idac=3.27e-07, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y',"Idac=1e-09, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=9e-08) Y',"Idac=9e-08, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y',"Idac=3.27e-07, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=1e-09) Y',"Idac=1e-09, LR mode")
    Graph = plot('servo_prelayout', 'Servo function - Pre layout, Vclamp=1.65', 'semilogx', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'Servo Magnitude', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    del Trace4, Trace4sel


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#SIZING OF THE PSEUDORESISTANCE ELEMENTS comparison    #################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if select:
    htmlPage("Pseudoresistor elements sizing")
    Trace1 = csv2traces('pseudoresistanceIdeal8BlocksModifiedW_version1.csv') #Ideal Vc
    Trace2 = csv2traces('pseudoresistanceIdeal8BlocksModifiedW_version2.csv')   #Ideal Ibias
    Trace3 = csv2traces('pseudoresistanceIdeal8Blocks.csv') #Ideal Vc
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['ExperimentCompounREq (Vc=0.493) Y', 'ExperimentCompounREq (Vc=0.138) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("ExperimentCompounREq (Vc=0.493) Y","Vc=0.493. W= 30u, L =500n")
            Trace1sel[key].label=Trace1sel[key].label.replace("ExperimentCompounREq (Vc=0.138) Y","Vc=0.138. W= 30u, L =500n")
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ['ExperimentCompounREq (Vc=0.398) Y', 'ExperimentCompounREq (Vc=0.002) Y']:
            print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("ExperimentCompounREq (Vc=0.398) Y","Vc=0.398. W= 30u, L =300n")
            Trace2sel[key].label=Trace2sel[key].label.replace("ExperimentCompounREq (Vc=0.002) Y","Vc=0.002. W= 30u, L =300n")
    Trace1sel.update(Trace2sel)
    keys=list(Trace3.keys())
    Trace3sel={}
    for key in keys:
        if key in ['ExperimentCompounREq (Vc=0.273) Y','ExperimentCompounREq (Vc=0.897) Y']:
            print(key)
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("ExperimentCompounREq (Vc=0.273) Y","Vc=0.273. W= 350n, L =500n (Selected)")
            Trace3sel[key].label=Trace3sel[key].label.replace("ExperimentCompounREq (Vc=0.897) Y","Vc=0.897. W= 350n, L =500n (Selected)")
    Trace1sel.update(Trace3sel)
    Graph = plot('WprLprSizing', 'Sizing of Wpr & Lpr. DC response.', 'semilogy', Trace1sel, xName = 'Vin', xUnits = 'V', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    del Trace1, Trace2, Trace3, Trace1sel, Trace2sel, Trace3sel

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#SIZING OF THE BIASING ELEMENTS   ######################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("Biasing elements sizing")
    Trace1 = Cadencecsv2traces('VbiasPseudoresistor.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ["/vbias (tryw=3.5e-07,tryl=5.175e-06) Y", "/vbias (tryw=1e-05,tryl=3.5e-07) Y"]:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("/vbias (tryw=3.5e-07,tryl=5.175e-06) Y","Wb1= 350n, Lb1=5u")
            Trace1sel[key].label=Trace1sel[key].label.replace("/vbias (tryw=1e-05,tryl=3.5e-07) Y","Wb2= 10u, Lb2=350n")
    Graph = plot('WbLbSizing', 'Sizing of Wb & Lb. DC response.', 'semilogx', Trace1sel, xName = 'Ic', xUnits = 'A', yName = 'Vc', yUnits = 'V', xLim = [1e-9,329e-9] , yLim = [], show = False)
    fig2html(Graph, 600)



########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#TRIODE VS CONVENTIONAL RESISTOR# ######################################################################################
select=True#############################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if select:
    htmlPage("Triode vs resistor")

    Trace = csv2traces('N_PMOS_triode_resistance.csv')
    Graph = plot('triodebehaviour', 'Equivalent resistance of the triode parallel resistor', 'semilogy', Trace, xName = '$V_{clamp}$', xUnits = 'V', yName = 'Resistance', yUnits = '$\Omega$', xLim = [1,2.3] , yLim = [1e5,1e12], show = False)
    fig2html(Graph, 600)

    Trace1 = Cadencecsv2traces('N_PMOS_triode_loopgain_mag_rph1.csv')
    Trace2 = Cadencecsv2traces('N_PMOS_triode_loopgain_mag_triode.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ["Loop Gain dB20 rnp1h (Vclamp=1.5,Rfb=5.023773e+08) Y"]:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain dB20 rnp1h (Vclamp=1.5,Rfb=5.023773e+08) Y","Resistor")
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ["Loop Gain dB20 triode (Vclamp=1.5,Rfb=5.023773e+08) Y"]:
            print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("Loop Gain dB20 triode (Vclamp=1.5,Rfb=5.023773e+08) Y","Triode")
    Trace1sel.update(Trace2sel)
    Graph = plot('triodevsresistorMag', 'Loop Gain Mag: rnp1h vs triode parallel resistors', 'semilogx', Trace1sel, xName = 'Frequency', xUnits = 'Hz', yName = 'Loop Gain Mag', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)

    Trace1 = Cadencecsv2traces('N_PMOS_triode_loopgain_phase_rph1.csv')
    Trace2 = Cadencecsv2traces('N_PMOS_triode_loopgain_phase_triode.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ["Loop Gain Phase rnp1h (Vclamp=1.5,Rfb=5.023773e+08) Y"]:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain Phase rnp1h (Vclamp=1.5,Rfb=5.023773e+08) Y","Resistor")
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ["Loop Gain Phase triode (Vclamp=1.5,Rfb=5.023773e+08) Y"]:
            print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("Loop Gain Phase triode (Vclamp=1.5,Rfb=5.023773e+08) Y","Triode")
    Trace1sel.update(Trace2sel)
    Graph = plot('triodevsresistorPhs', 'Loop Gain Phase: rnp1h vs triode parallel resistors', 'semilogx', Trace1sel, xName = 'Frequency', xUnits = 'Hz', yName = 'Loop Gain Phase', yUnits = '$deg$', xLim = [1e-2,10e10] , yLim = [-200,200], show = False)
    fig2html(Graph, 600)
"""
