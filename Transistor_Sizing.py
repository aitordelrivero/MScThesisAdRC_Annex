#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from SLiCAP import *
from sympy.plotting import plot3d
from scipy.optimize import fminbound
from math import ceil, prod, log

def plotSweepEdited(fileName, title, results, sweepStart, sweepStop, sweepNum, sweepVar = 'auto', sweepScale = '', xVar = 'auto', xScale = '', xUnits = '', axisType = 'auto', funcType = 'auto', yVar = 'auto', yScale = '', yUnits = '', noiseSources = None, show = False):
    """
    THE FUNCTION IS ALMOST EQUAL TO THE ORIGINAL IN SLICAP.

    Plots a function by sweeping one variable and optionally stepping another.
    The function to be plotted depends on the arguments 'yVar' and 'funcType':
    - If funcType == 'params', the variable 'yVar' must be the name of a circuit
      parameter.
    - If funcType == 'auto', the default function that will be plotted depends
      on the data type of the instruction:
      - data type == 'noise': funcType = 'onoise'
      - data type == 'laplace', 'numer' or 'denom': funcType = 'mag'
      - data type == 'time', 'impulse' or 'step': funcType = 'time'
    The variable plotted along the x-axis defaults to the sweep variable. However,
    for multivariate functions obtained with data type 'params', the x variable
    can be choosen from all circuit parameters.
    - If sweepVar == 'auto', the sweep variable will be determined from the data type:
      - data type == 'noise', 'laplace', 'numer' or 'denom': sweepVar = ini.frequency
        for data types 'laplace', 'numer' or 'denom' the laplace variable will
        be replaced with sympy.i*ini.frequency or with 2*sympy.pi*sympy.i*ini.frequency
        before sweeping, when ini.Hz == False, or ini.Hz== True, respectively.
      - dataType == 'time', 'impulse' or 'step': sweepVar = sympy.Symbol('t')
    The type of axis can be 'lin', 'log', 'semilogx', 'semilogy' or 'polar'.
    :param fileName: Name of the file for saving it to disk.
    :type fileName: str
    :param title: Title of the figure.
    :type title: str
    :param results: Results of the execution of an instruction, or a list with
                    SLiCAPprotos.allResults objects.
    :type results: list, SLiCAPprotos.allResults
    :param sweepStart: Start value of the sweep parameter
    :type sweepStart: float, int, str
    :param sweepStop: Stop value of the sweep parameter
    :type sweepStop: float, int, str
    :param sweepNum: Number of points of the sweep parameter
    :type sweepNum: int
    :param sweepVar: Name of the sweep variable
    :type sweepVar: sympy.Symbol, str
    :param sweepScale: Scale factor of the sweep variable. Both the start and
                       the stop value will be multiplied with a factor that
                       corresponds with this scale factor.
    :type sweepScale: str
    :param xVar: Name of the variable to be plotted along the x axis
    :type xvar: str, sympy.Symbol
    :param xScale: Scale factor of the x axis variable.
    :type xScale: str
    :param xUnits: Units of the x axis variable.
    :type xUnits: str
    :param axisType: Type of axis: 'lin', 'log', 'semilogx', 'semilogy' or 'polar'.
    :type axisType: str
    :param funcType: Type of function can be: 'mag', 'dBmag', 'phase', 'delay',
                     'time', 'onoise', 'inoise' or 'param'.
    :type funcType: str
    :param yScale: Scale factor of the y axis variable.
    :type yScale: str
    :param yUnits: Units of the y axis variable.
    :type yUnits: str
    :param noiseSources: Noise sources of which the contribution to the detector-
                         referred noise (funcType = 'onoise') or the source-
                         referred noise (funcType = 'inoise') should be plotted.
                         Can be 'all', a list with names of noise sources or an
                         ID of a noise source.
    :type noiseSources: list, str
    :param show: If 'True' the plot will be shown in the workspace.
    :type show: bool
    :return: fig
    :rtype: SLiCAPplots.figure
    """
    ini.defaultColors=['#c3312f','#00b7d3','#99d28c','#f1be3e','#007188','k','#eb7245','#0066a2','#00a390','#82d7c6']
    ini.plotFontSize=18
    plotDataTypes = ['laplace', 'numer', 'denom', 'noise', 'step', 'impulse', 'time', 'params', None]
    funcTypes  = ['mag', 'dBmag', 'phase', 'delay', 'time', 'onoise', 'inoise', 'param']
    axisTypes  = ['lin', 'log', 'semilogx', 'semilogy', 'polar']
    freqTypes  = ['laplace', 'numer', 'denom', 'noise']
    timeTypes  = ['time', 'impulse', 'step']
    fig = figure(fileName)
    fig.show = show
    ax = axis(title)
    if type(results) != list:
        results = [results]
    colNum = 0
    numColors = len(ini.defaultColors)
    # first results defines the axis type and labels
    result = results[0]
    if result.dataType not in plotDataTypes:
        print("Error: cannot plot dataType '{0}' with 'plotSweep()'.".format(result.dataType))
        return fig
    if funcType == 'auto':
        if result.dataType == 'noise':
            funcType = 'onoise'
        elif result.dataType in freqTypes:
            funcType = 'mag'
        elif result.dataType in timeTypes:
            funcType = 'time'
    elif funcType == 'param':
        if sweepVar == 'auto':
            print("Error: undefined sweep variable.")
            return fig
        if yVar == 'auto':
            print("Error: missing parameter to be plotted.")
            return fig
        if xVar == 'auto':
            xVar = sweepVar
            xScale = sweepScale
    elif funcType not in funcTypes:
        print("Error: unknown funcType: '{0}'.".format(funcType))
        return fig
    if axisType == 'auto':
        if funcType == 'param':
            axisType = 'lin'
        elif funcType == 'mag' or result.dataType == 'noise':
            axisType = 'log'
        elif funcType == 'dBmag' or funcType == 'phase' or funcType == 'delay':
            axisType = 'semilogx'
        elif funcType == 'time':
            axisType = 'lin'
    elif axisType not in axisTypes:
        print("Error: unknown axisType: '{0}'.".format(axisType))
        return fig
    if axisType == 'lin':
        ax.xScale = 'lin'
        ax.yScale = 'lin'
    elif axisType == 'log':
        ax.xScale = 'log'
        ax.yScale = 'log'
    elif axisType == 'semilogx':
        ax.xScale = 'log'
        ax.yScale = 'lin'
    elif axisType == 'semilogy':
        ax.xScale = 'lin'
        ax.yScale = 'log'
    elif axisType == 'polar':
        ax.polar = True
        ax.yScale = 'lin'
    if not ax.polar:
        if funcType == 'param' and xVar != sweepVar:
            ax.xScaleFactor = xScale
        else:
            ax.xScaleFactor = sweepScale
    ax.yScaleFactor = yScale
    ax.traces = []
    if funcType == 'param':
        #####################################################
        if str(sp.latex(sp.Symbol(xVar)))=='Lcm':
            ax.xLabel = '$L_{M2,M3}$[' + xScale + xUnits + ']'
        elif str(sp.latex(sp.Symbol(xVar)))=='Lbcm':
            ax.xLabel = '$L_{M6}$[' + xScale + xUnits + ']'
        elif str(sp.latex(sp.Symbol(xVar)))=='Lb':
            ax.xLabel = '$L_{M5}$[' + xScale + xUnits + ']'
        elif str(sp.latex(sp.Symbol(xVar)))=='Lout':
            ax.xLabel = '$L_{M4}$[' + xScale + xUnits + ']'
        else:
            ax.xLabel = '$' + sp.latex(sp.Symbol(xVar)) + '$ [' + xScale + xUnits + ']'

        if str(sp.latex(sp.Symbol(yVar)))=='f_{T XU3}':
            ax.yLabel = '$f_{T(M2,M3)}$[' + yScale + yUnits + ']'
        elif str(sp.latex(sp.Symbol(yVar)))=='MaxIN':
            ax.yLabel = '$V_{in,max}$[' + yScale + yUnits + ']'
        elif str(sp.latex(sp.Symbol(yVar)))=='VIN_{min}':
            ax.yLabel = '$V_{in,min}$[' + yScale + yUnits + ']'
        elif str(sp.latex(sp.Symbol(yVar)))=='Vout_{min}':
            ax.yLabel = '$V_{out,min}$[' + yScale + yUnits + ']'
        elif str(sp.latex(sp.Symbol(yVar)))=='Vout_{max}':
            ax.yLabel = '$V_{out,max}$[' + yScale + yUnits + ']'
        else:
            ax.yLabel = '$' + sp.latex(sp.Symbol(yVar)) + '$[' + yScale + yUnits + ']'
            #####################################################
    elif result.dataType in freqTypes:
        if result.dataType == 'noise':
            if funcType == 'onoise':
                yUnits = result.detUnits
            if funcType == 'inoise':
                yUnits = result.srcUnits
            ax.xLabel = 'frequency [' + sweepScale + 'Hz]'
        elif ini.Hz == True:
            ax.xLabel = 'frequency [' + sweepScale + 'Hz]'
        else:
            ax.xLabel = 'frequency [' + sweepScale + 'rad/s]'
    # For time plots we use time along the x-axis
    elif funcType in timeTypes:
        ax.xLabel = 'time [' + sweepScale + 's]'
        yUnits = result.detUnits
    # Create the y-label for other than parameter plots
    if funcType == 'mag':
        ax.yLabel = 'magnitude [' + yScale + yUnits + ']'
    elif funcType == 'dBmag':
        ax.yLabel = 'magnitude [' + yScale + 'dB]'
    elif funcType == 'phase':
        if ini.Hz == True:
            ax.yLabel = 'phase [' + yScale + 'deg]'
        else:
            ax.yLabel = 'phase [' + yScale + 'rad]'
    elif funcType == 'delay':
        ax.yLabel = 'group delay [' + yScale + 's]'
    elif funcType == 'time':
        ax.yLabel = '[' + yScale + yUnits + ']'
    elif funcType == 'onoise':
        ax.yLabel = 'spectral density [$' + yScale + yUnits +'^2/Hz$]'
    elif funcType == 'inoise':
        ax.yLabel = 'spectral density [$' + yScale + yUnits +'^2/Hz$]'
    # Create the sweep, lin or log depending on the x-axis type
    try:
        xScaleFactor = 10**int(SCALEFACTORS[sweepScale])
    except:
        xScaleFactor = 1.
    if ax.xScale == 'log' or ax.xScale == 'semilogx' or ax.polar == True:
        x = np.geomspace(checkNumber(sweepStart)*xScaleFactor, checkNumber(sweepStop)*xScaleFactor, checkNumber(sweepNum))
    elif ax.xScale == 'lin' or ax.xScale == 'semilogy':
        x = np.linspace(checkNumber(sweepStart)*xScaleFactor, checkNumber(sweepStop)*xScaleFactor, checkNumber(sweepNum))
    if funcType == 'param':
        xData, yData = stepParams(result, xVar, yVar, sweepVar, x)
        if type(xData) == dict:
            keys = sorted(list(xData.keys()))
            for i in range(len(keys)):
                newTrace = trace([xData[keys[i]], yData[keys[i]]])
                #####################################################
                if str(result.stepVar)=='Wcm':
                    newTrace.label = '$W_{M2,M3}$ = %8.1e'%(result.stepList[i])
                elif str(result.stepVar)=='Wbcm':
                    newTrace.label = '$W_{M6}$ = %8.1e'%(result.stepList[i])
                elif str(result.stepVar)=='Wb':
                    newTrace.label = '$W_{M5}$ = %8.1e'%(result.stepList[i])
                elif str(result.stepVar)=='Wout':
                    newTrace.label = '$W_{M4}$ = %8.1e'%(result.stepList[i])
                else:
                    newTrace.label = '$%s$ = %8.1e'%(sp.latex(result.stepVar), result.stepList[i])
                #####################################################
                newTrace.color = ini.defaultColors[colNum % numColors]
                colNum += 1
                print(newTrace.label)
                newTrace.label.replace("loopgain, ","")
                ax.traces.append(newTrace)
        else:
            newTrace = trace([xData, yData])
            newTrace.label = '$' + sp.latex(sp.Symbol(yVar)) + '$'
            newTrace.color = ini.defaultColors[colNum % numColors]
            ax.traces.append(newTrace)
            colNum += 1
    else:
        for result in results:
            if not result.step:
                if result.dataType == 'numer':
                    yData = result.numer
                    yLabel = 'numer: '
                elif result.dataType == 'denom':
                    yData = result.denom
                    yLabel = 'denom: '
                elif result.dataType == 'laplace':
                    yData = result.laplace
                    yLabel = ''
                elif result.dataType == 'time':
                    yData = result.time
                    yLabel = '$' + sp.latex(sp.Symbol(result.detLabel)) + '$'
                elif result.dataType == 'step':
                    yData = result.stepResp
                    yLabel = '$' + sp.latex(sp.Symbol(result.detLabel)) + '$'
                elif result.dataType == 'impulse':
                    yData = result.impulse
                    yLabel = '$' + sp.latex(sp.Symbol(result.detLabel)) + '$'
                if funcType == 'mag':
                    if ax.polar:
                        radius = magFunc_f(yData, x)
                        angle = phaseFunc_f(yData, x)
                        if ini.Hz:
                            angle = angle/180*np.pi
                        newTrace = trace([angle, radius])
                    else:
                        newTrace = trace([x, magFunc_f(yData, x)])
                elif funcType == 'dBmag':
                    if ax.polar:
                        radius = dBmagFunc_f(yData, x)
                        angle = phaseFunc_f(yData, x)
                        if ini.Hz:
                            angle = angle/180*np.pi
                        newTrace = trace([angle, radius])
                    else:
                        newTrace = trace([x, dBmagFunc_f(yData, x)])
                elif funcType == 'phase':
                    if not ax.polar:
                        newTrace = trace([x, phaseFunc_f(yData, x)])
                elif funcType == 'delay':
                    if not ax.polar:
                        newTrace = trace([x, delayFunc_f(yData, x)])
                elif funcType == 'time':
                    if not ax.polar:
                        if sp.Symbol('t') in list(yData.atoms(sp.Symbol)):
                            func = sp.lambdify(sp.Symbol('t'), yData)
                            y = np.real(func(x))
                        else:
                            y = [yData for i in range(len(x))]
                        newTrace = trace([x, y])
                if result.dataType != 'noise':
                    if result.gainType == 'vi':
                        newTrace.label = result.detLabel
                    else:
                        newTrace.label = result.gainType
                    try:
                        newTrace.color = ini.gainColors[result.gainType]
                    except:
                        newTrace.color = ini.defaultColors[colNum % numColors]
                        colNum += 1
                    if result.gainType == 'vi':
                        yLabel += '$' + sp.latex(sp.Symbol(result.detLabel)) + '$'
                    else:
                        yLabel += result.gainType
                    try:
                        ax.traces.append(newTrace)
                    except:
                        pass
                else:
                    keys = list(result.onoiseTerms.keys())
                    if noiseSources == None:
                        if funcType == 'onoise':
                            yData = result.onoise
                        elif funcType == 'inoise':
                            yData = result.inoise
                        func = sp.lambdify(ini.frequency, yData)
                        # y = [func(x[j]) for j in range(len(x))]
                        if ini.frequency in list(yData.atoms(sp.Symbol)):
                            y = func(x)
                        else:
                            y = [yData for i in range(len(x))]
                        newTrace = trace([x, y])
                        newTrace.label = funcType
                        ax.traces.append(newTrace)
                    elif noiseSources == 'all':
                        for srcName in keys:
                            if funcType == 'onoise':
                                yData = result.onoiseTerms[srcName]
                            elif funcType == 'inoise':
                                yData = result.inoiseTerms[srcName]
                            func = sp.lambdify(ini.frequency, yData)
                            # y = [func(x[j]) for j in range(len(x))]
                            if ini.frequency in list(yData.atoms(sp.Symbol)):
                                y = func(x)
                            else:
                                y = [yData for i in range(len(x))]
                            noiseTrace = trace([x, y])
                            noiseTrace.color = ini.defaultColors[colNum % numColors]
                            noiseTrace.label = funcType + ': ' + srcName
                            ax.traces.append(noiseTrace)
                            colNum += 1
                    elif noiseSources in keys:
                        if funcType == 'onoise':
                            yData = result.onoiseTerms[noiseSources]
                        elif funcType == 'inoise':
                            yData = result.inoiseTerms[noiseSources]
                        func = sp.lambdify(ini.frequency, yData)
                        # y = [func(x[j]) for j in range(len(x))]
                        if ini.frequency in list(yData.atoms(sp.Symbol)):
                            y = func(x)
                        else:
                            y = [yData for i in range(len(x))]
                        noiseTrace = trace([x, y])
                        noiseTrace.color = ini.defaultColors[colNum % numColors]
                        noiseTrace.label = funcType + ': ' + sources
                        ax.traces.append(noiseTrace)
                        colNum += 1
                    elif type(noiseSources) == list:
                        for srcName in noiseSources:
                            if srcName in keys:
                                if funcType == 'onoise':
                                    yData = result.onoiseTerms[srcName]
                                elif funcType == 'inoise':
                                    yData = result.inoiseTerms[srcName]
                                func = sp.lambdify(ini.frequency, yData)
                                # y = [func(x[j]) for j in range(len(x))]
                                if ini.frequency in list(yData.atoms(sp.Symbol)):
                                    y = func(x)
                                else:
                                    y = [yData for i in range(len(x))]
                                noiseTrace = trace([x, y])
                                noiseTrace.color = ini.defaultColors[colNum % numColors]
                                noiseTrace.label = funcType + ': ' + srcName
                                ax.traces.append(noiseTrace)
                                colNum += 1
                    else:
                        print("Error: cannot understand 'sources={0}'.".format(str(sources)))
                        return fig
            else:
                if result.stepMethod != 'array':
                    stepNum = len(result.stepList)
                else:
                    stepNum = len(result.stepArray[0])
                for i in range(stepNum):
                    if result.dataType == 'numer':
                        yData = result.numer[i]
                        yLabel = 'numer: '
                    elif result.dataType == 'denom':
                        yData = result.denom[i]
                        yLabel = 'denom: '
                    elif result.dataType == 'laplace':
                        yData = result.laplace[i]
                        yLabel = ''
                    elif result.dataType == 'time':
                        yData = result.time[i]
                    elif result.dataType == 'step':
                        yData = result.stepResp[i]
                    elif result.dataType == 'impulse':
                        yData = result.impulse[i]
                    elif result.dataType == 'noise':
                        if funcType == 'onoise':
                            yData = result.onoise[i]
                        elif funcType == 'inoise':
                            yData = result.inoise[i]
                    if result.gainType == 'vi':
                        if result.dataType == 'noise':
                            yLabel = funcType
                        else:
                            yLabel += '$' + sp.latex(sp.Symbol(result.detLabel)) + '$'
                    else:
                        yLabel = result.gainType
                    if result.stepMethod == 'array':
                        yLabel += ', run: %s'%(i+1)
                    else:
                        yLabel += ', %s = %8.1e'%(result.stepVar, result.stepList[i])
                    if funcType == 'mag':
                        if ax.polar:
                            radius = magFunc_f(yData, x)
                            angle = phaseFunc_f(yData, x)
                            if ini.Hz:
                                angle = angle/180*np.pi
                            newTrace = trace([angle, radius])
                        else:
                            newTrace = trace([x, magFunc_f(yData, x)])
                    elif funcType == 'dBmag':
                        if ax.polar:
                            radius = dBmagFunc_f(yData, x)
                            angle = phaseFunc_f(yData, x)
                            if ini.Hz:
                                angle = angle/180*np.pi
                            newTrace = trace([angle, radius])
                        else:
                            newTrace = trace([x, dBmagFunc_f(yData, x)])
                    elif funcType == 'phase':
                        if not ax.polar:
                            newTrace = trace([x, phaseFunc_f(yData, x)])
                    elif funcType == 'delay':
                        if not ax.polar:
                            newTrace = trace([x, delayFunc_f(yData, x)])
                    elif funcType == 'time':
                        if not ax.polar:
                            if sp.Symbol('t') in list(yData.atoms(sp.Symbol)):
                                func = sp.lambdify(sp.Symbol('t'), yData)
                                y = np.real(func(x))
                            else:
                                y = [yData for i in range(len(x))]
                            newTrace = trace([x, y])
                    elif funcType == 'onoise' or funcType == 'inoise':
                        if not ax.polar:
                            func = sp.lambdify(ini.frequency, yData)
                            #y = [func(x[j]) for j in range(len(x))]
                            if ini.frequency in list(yData.atoms(sp.Symbol)):
                                y = func(x)
                            else:
                                y = [yData for i in range(len(x))]
                            newTrace = trace([x, y])
                    newTrace.color = ini.defaultColors[colNum % numColors]
                    colNum += 1
                    newTrace.label = yLabel
                    try:
                        alternativeTag=newTrace.label.replace("loopgain, "," ")
                        newTrace.label=alternativeTag
                        ax.traces.append(newTrace)
                    except:
                        pass

    fig.axes = [[ax]]
    fig.plot()
    ini.defaultColors      = ['k','b','g','c','m','y','k']
    ini.plotFontSize=12
    return fig


t1 = time()
prj = initProject('MSc Thesis')
print("Project inititated")
ini.MaximaTimeOut = 600         #Some extra time to allow the computer to calculate integrals.


# Transistor sizing (Instance i1)
Wopt=30e-6
Iopt=2e-6
L = 2e-6
Wout=5e-6
Lout=350e-9
IDout=5e-6
Lcm=5e-6
Wcm=7.5e-6
Wb=5e-6
Lb=1.5e-6
Wbcm=8e-6
Lbcm=1.5e-6

CpsDefaultBIG='5p' #12p
CpsDefaultSMALL='5p' #3p
RpsDefault=200000 #None for 1/gm4
Rph1Default=4000

CompensationArea=26.5e-6**2+(4+1)*350e-9*10e-6+1e-6*10e-6 #Cps+RpsPMOS+RpsNMOS+Rph1 areas

CparX=10.5e-15
CparC=25.03e-15
CparY=29.73e-15


#Model parameters and specifications
IG = 0                          #Gate current is neglectible for optimization
Ce = 5e-12                      #Neuron and electrode impedances values are calculated as explained in the introduction (Wiki).
Ree = 1500e9
Rs = 100e6
Cjm = 0.1e-12
Rjm = 1e9
Rm = 100e6
Cm = 100e-12
f_min = 10e-4                     #Hz (ideally DC, but some margin to avoid results involving infinity)
f_max = 10e3                      #Hz (maximum frequency of the input signal)
max_gain = 20e9                   #Maxiumum gain of the TIA as explained in the amplifier type page (Wiki)
min_gain = 2e6                    #Minimum gain of the TIA as explained in the amplifier type page (Wiki)
m_gain = 20e7                     #An intermediate gain for checking purposes.

Cadc=5e-12
Cgm=0
Cload=Cadc+Cgm

Vdd=3.3

#ADC specs
res=10                          #-bit
rangeADC=2                         #V
f_sampling=20e3                 #sampling frequency of the ADC

#Execution parameters:
ngtbc=3 #Number of gains to be checked
GteChnNoiseFactor=9


#List of all circuits that will be checked
circuit="FinalCircuit" #"AllParasitics"
Name=circuit
test=False

makeNetlist(Name+'.asc',Name)       #create netlist
i1=instruction()                    #create instruction object
i1.setCircuit(Name+'.cir')          #load circuit
print("-------------------")
print("CIRCUIT SET: "+Name)


#Loading model parameters into the .cir
#######???????????????????????????????
#COPY HERE PROCESS PARAMETERS
#######???????????????????????????????


if (Name=="CIRCUIT1" or Name=="SingleStage" or Name=="TopBiasing1" or Name=="BottomBiasing1" or Name=="TwoStages" or Name=="FreqComp_Both_NoBias" or Name=="Biasing2" or Name=="FinalCircuit" or Name=="ParasiticCM" or Name=="ParasiticNET108" or Name=="ParasiticOUT1" or Name=="AllParasitics"):
    reference='Gm_M1_XU1'
    i1.defPar('c_dg_XU1', 0)
    #i2.defPar('c_dg_XU1', 0)
    print('loop gain reference assigned, and local loop in main signal path deleted')
if (Name=="TopBiasing1" or Name=="BottomBiasing1" or Name=="TwoStages" or Name=="Biasing2" or Name=="FinalCircuit" or Name=="FreqComp_Both_NoBias" or Name=="ParasiticCM" or Name=="ParasiticNET108" or Name=="ParasiticOUT1" or Name=="AllParasitics"):
    i1.defPar('c_dg_XU3', 0)
    i1.defPar('c_dg_XU2', 0)
    # i2.defPar('c_dg_XU3', 0)
    print('local loop in secondary signal path deleted')
else:
    print('loop gain reference was not assigned')
try:
    i1.defPar('W', Wopt)
    i1.defPar('ID', Iopt)
    i1.defPar('Lin', L)
    print('single-stage parameters assigned')
except:
    print('single-stage parameters not needed')
try:
    i1.defPar('IDoutPMOS', -IDout)
    i1.defPar('Wout', Wout)
    i1.defPar('Lout', Lout)
    i1.defPar('Vdd', Vdd)
    i1.defPar('range_ADC', rangeADC)
    print('output parameters assigned')
except:
    print('two-stage ouput parameters not needed')
try:
    i1.defPar('IDas', Iopt*2)
    i1.defPar('Was', Wopt*2)
    i1.defPar('L', L)
    print('input balanced CS parameters assigned')
except:
    print('input antiseries not needed')
try:
    i1.defPar('Wcm', Wcm)
    i1.defPar('Lcm', Lcm)
    i1.defPar('IDasPMOS', -Iopt*2)
    print('current mirror parameters assigned')
except:
    print('no current mirrors')
try:
    i1.defPar('Wb', Wb)
    i1.defPar('Lb', Lb)
    i1.defPar('IDout', IDout)
    print('output bias equivalent parameters assigned')
except:
    print('no output bias equivalent')
try:
    i1.defPar('Wbcm', Wbcm)
    i1.defPar('Lbcm', Lbcm)
    i1.defPar('IDas', Iopt*2)
    print('output bias equivalent parameters assigned')
except:
    print('no input bias equivalent')
try:
    i1.defPar('C_par_CM', CparX)
    print('parasitic C_par_CM set')
except:
    print('parasitic C_par_CM not modeled')
try:
    i1.defPar('C_par_OUT1', CparC)
    print('parasitic C_par_OUT1 set')
except:
    print('parasitic C_par_OUT1 not modeled')
try:
    i1.defPar('C_par_NET018', C_par_NET018)
    print('parasitic C_par_NET018 set')
except:
    print('parasitic C_par_NET018 not modeled')


#t_ox, Vth, n, u0, E_crit and theta, for both NMOS and PMOS is defined locally since it is confidential data.
i1.defPar('IG', IG)
i1.defPar('Ce', Ce)
i1.defPar('Ree', Ree)
i1.defPar('Rs', Rs)
i1.defPar('Cload', Cload)
print("Parameters assigned")


gains=np.geomspace(min_gain, max_gain, num=ngtbc, endpoint=True).tolist()#generate list of gains to be checked

#Print all circuit data
htmlPage('Circuit data')
head2html('Circuit diagram')
img2html(Name + '.png', 1200)
print("Circuit data printed")
netlist2html(Name + '.cir')

for gain in gains:
    try:
        gm4=i1.getParValue('g_m_XU4')
    except:
        print("no gm4")
    try:
        if RpsDefault==None:
            i1.defPar('Rps', 1/gm4)
        else:
            i1.defPar('Rps', RpsDefault)
    except:
        print('no Rps')
    if gain<200e6:
        try:
            i1.defPar('Cps', CpsDefaultBIG)
        except:
            print('no Cps')
    else:
        try:
            i1.defPar('Cps', CpsDefaultSMALL)
        except:
            print('no Cps')
    try:
        i1.defPar('Rph1', Rph1Default)
    except:
        print('no Rph1')

    i1.defPar('gain', gain)
    head2html('Gain='+str(format(gain,'.1E')))
    head2html('Element Data')
    elementData2html(i1.circuit)
    head2html('Circuit Parameters')
    params2html(i1.circuit)
    i1.setSimType('numeric') #numeric simulation to substitute all variables with their values
    i1.setSource('I1')       #the source is the membrane current source
    i1.setDetector('V_out')  #the detector is the output (input of the ADC)
    i1.setLGref(reference) #reference source in the asymptotic model
    i1.setGainType('loopgain')
    i1.setDataType('pz')
    loopgainpz=i1.execute()
    pz2html(loopgainpz) #DEBUG: print pole zero results
    i1.setDataType('poles')
    loopgainpoles=i1.execute()
    pz2html(loopgainpoles) #DEBUG: print pole zero results
    i1.setDataType('zeros')
    loopgainzeros=i1.execute()
    pz2html(loopgainzeros) #DEBUG: print pole zero results
    i1.setGainType('servo')
    i1.setDataType('pz')
    servopz=i1.execute()
    pz2html(servopz) #DEBUG: print pole zero results


#     #ASYMPTOTIC GAIN AND POLE AND ZERO ANALYSIS
f_check=10e6
i1.setSimType('numeric') #numeric simulation to substitute all variables with their values
i1.setSource('I1')       #the source is the membrane current source
i1.setDetector('V_out')  #the detector is the output (input of the ADC)
i1.setLGref(reference) #reference source in the asymptotic model
i1.setStepVar('gain')
i1.setStepMethod('log')
i1.setStepStart(min_gain)
i1.setStepStop(max_gain)
i1.setStepNum(ngtbc)
i1.stepOn()

i1.setGainType('loopgain')
i1.setDataType('laplace')
loopgain = i1.execute()
i1.setGainType('servo')
servo=i1.execute()
i1.setGainType('gain')
gain=i1.execute()
i1.setGainType('asymptotic')
asymptotic=i1.execute()
i1.setGainType('direct')
direct=i1.execute()


#print(loopgain.laplace[0])
ini.plotFontSize=18
head2html('Bode plot. Loopgain.')
result = [loopgain]
figdBmag = plotSweepEdited('dBmagL'+Name+'SizingScript', 'Loop gain - dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
figPhase = plotSweepEdited('phaseL'+Name+'SizingScript', 'Loop gain - Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
fig2html(figdBmag, 800)
fig2html(figPhase, 800)
print("Loopgain bode obtained")

#print(servo.laplace[0])
head2html('Bode plot. Servo.')
result = [servo]
figdBmag = plotSweep('figdBmagS'+Name+'SizingScript', 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
figPhase = plotSweep('phaseS'+Name+'SizingScript', 'Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
fig2html(figdBmag, 800)
fig2html(figPhase, 800)
print("Servofunction bode obtained")

#print(gain.laplace[0])
head2html('Bode plot. Gain.')
result = [gain]
figdBmag = plotSweep('figdBmagG'+Name+'SizingScript', 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
figPhase = plotSweep('phaseG'+Name+'SizingScript', 'Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
fig2html(figdBmag, 800)
fig2html(figPhase, 800)
print("Gain bode obtained")

#print(asymptotic.laplace[0])
head2html('Bode plot. Asymptotic.')
result = [asymptotic]
figdBmag = plotSweep('figdBmagA'+Name+'SizingScript', 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
figPhase = plotSweep('phaseA'+Name+'SizingScript', 'Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
fig2html(figdBmag, 800)
fig2html(figPhase, 800)
print("Asymptotic gain bode obtained")


#print(direct.laplace[0])
head2html('Bode plot. Direct.')
result = [direct]
figdBmag = plotSweep('figdBmagD'+Name+'SizingScript', 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
figPhase = plotSweep('phaseD'+Name+'SizingScript', 'Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
fig2html(figdBmag, 800)
fig2html(figPhase, 800)
print("Direct gain bode obtained")

i1.stepOff()


if Name=="TopBiasing1":
    ini.defaultColors=['#c3312f','#00b7d3','#99d28c','#f1be3e','#007188','k','#eb7245','#0066a2','#00a390','#82d7c6']
    #Noise performance
    i1.defPar('NEF', '1+(g_m_XU3/(g_m_XU1))') #"Noise excess factor"
    #Gain performance
    i1.defPar('Av_stg1', '(g_m_XU1)*(1/((g_o_XU1)+g_o_XU3))*(2*g_m_XU3*(1/g_o_XU3)+1)/(2*(g_m_XU3*(1/g_o_XU3)+1))') #i expect almost twice the next (always between 1.5 and 2 times)
    i1.defPar('Av_stg1_woCM', '((g_m_XU1)/2)*(1/((g_o_XU1)+g_o_XU3))') #considering that the biasing transistor is like the right current mirror transistor (i expect less than 1400, let's say X)
    i1.defPar('Av_stg1_woCMideal', '(g_m_XU1)/(g_o_XU1)') #considering circuit 2 (expected gain: 1400 approx)
    i1.defPar('GEFcm', '(Av_stg1)/Av_stg1_woCMideal') #"Gain excess factor"
    #Noise-Gain performance
    i1.defPar('GONcm', 'GEFcm/NEF') #"Gain excess factor over Noise excess factor"
    #Area performance
    i1.defPar('AreaCM', '2*Wcm*Lcm')
    i1.defPar('AreaBCS', '2*Was*L')
    i1.defPar('AREAfactor', '(AreaCM+AreaBCS)/(AreaBCS)')
    #Noise-gain-area performance
    i1.defPar('GONAcm', 'GEFcm/(NEF*AREAfactor)')
    #Input sWoptg performance
    i1.defPar('Vsat_XU1','(2*U_T*sqrt(IC_XU1+0.25)+3*U_T)/2')
    i1.defPar('Vsat_XU2','2*U_T*sqrt(IC_XU2+0.25)+3*U_T')
    i1.defPar('Vsat_XU3','2*U_T*sqrt(IC_XU3+0.25)+3*U_T')
    i1.defPar('MaxIN_woCM', 'Vdd-Vsat_XU3-Vsat_XU1+V_GS_XU1')
    i1.defPar('MaxIN', 'Vdd+V_GS_XU3-Vsat_XU1+V_GS_XU1')
    i1.defPar('MaxINratio', 'MaxIN/MaxIN_woCM')
    #PSRR
    i1.defPar('PSRR','g_m_XU1/(g_o_XU1+g_o_XU3)')

    params2html(i1.circuit)

    i1.setDataType('params')
    i1.stepOn()
    i1.setStepVar('Wcm')
    i1.setStepStart(0.5e-6)
    i1.setStepStop(60e-6)
    i1.setStepNum(10)
    i1.setStepMethod('log')
    result = i1.execute()

    head2html('Adjusting Wcm based on ft requirement')
    Sweep = plotSweepEdited('WCurrentMirror1', '$Speed: f_T$', result, 350e-9, 10e-6, 50, axisType = 'semilogy', sweepVar= 'Lcm', xUnits = 'm', yVar = 'f_T_XU3', yUnits = 'Hz', funcType = 'param', show = False)
    fig2html(Sweep, 600)

    head2html('Adjusting Wcm based on noise requirement')
    Sweep2 = plotSweepEdited('WCurrentMirror2', 'Noise: NEF', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lcm', xUnits = 'm', yVar = 'NEF', yUnits = 'Hz', funcType = 'param', show = False)
    fig2html(Sweep2, 600)

    head2html('Adjusting Wcm based on gain requirement')
    Sweep3 = plotSweepEdited('WCurrentMirror3', '$GEF(L_{CurMir})$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lcm', xUnits = 'm', yVar = 'GEFcm', yUnits = '', funcType = 'param', show = False)
    fig2html(Sweep3, 600)
    text2html('The GOF is less than 1 because it is compared to the ideal situation in which the biasing is a current source with infinite output impedance. The GOF is around 2 if we consider that source.')

    head2html('Trade-off: GEF over NEF')
    Sweep4 = plotSweepEdited('WCurrentMirror4', '$GON(L_{CurMir})$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lcm', xUnits = 'm', yVar = 'GONcm', yUnits = '', funcType = 'param', show = False)
    fig2html(Sweep4, 600)


    head2html('Trade-off: GEF over NEF and area')
    Sweep5 = plotSweepEdited('WCurrentMirror5', '$GONA(L_{CurMir})$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lcm', xUnits = 'm', yVar = 'GONAcm', yUnits = '', funcType = 'param', show = False)
    fig2html(Sweep5, 600)


    head2html('PSRR')
    Sweep8 = plotSweepEdited('WCurrentMirror8', '$Isolation: PSRR$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lcm', xUnits = 'm', yVar = 'PSRR', yUnits = '', funcType = 'param', show = False)
    fig2html(Sweep8, 600)


    head2html('MaxIN')
    Sweep9 = plotSweepEdited('WCurrentMirror9', '$Swing: Vin_{max}$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lcm', xUnits = 'm', yVar = 'MaxIN', yUnits = 'V', funcType = 'param', show = False)
    fig2html(Sweep9, 600)
    text2html('Target: '+str(Vdd-rangeADC/2))

    head2html('Adjusting Wcm based on sWoptg requirement')
    Sweep6 = plotSweepEdited('WCurrentMirror6', '$MaxINratio(L_{CurMir})$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lcm', xUnits = 'm', yVar = 'MaxINratio', yUnits = '', funcType = 'param', show = False)
    fig2html(Sweep6, 600)

    i1.stepOff()

    ini.defaultColors      = ['k','b','g','c','m','y','k']


if Name=="BottomBiasing1":
    head2html('Trade-offs')
    #CMRR
    i1.defPar('CMRR', '(1+2*(g_m_XU1)*(1/g_o_XU6))*g_m_XU3*(1/(g_o_XU1+g_o_XU3))')
    #SWoptg
    i1.defPar('Vsat_XU1','(2*U_T*sqrt(IC_XU1+0.25)+3*U_T)/2')
    i1.defPar('Vsat_XU2','2*U_T*sqrt(IC_XU2+0.25)+3*U_T') #Binkley's appraoch
    i1.defPar('Vsat_XU3','2*U_T*sqrt(IC_XU3+0.25)+3*U_T')
    i1.defPar('Vsat_XU6','2*U_T*sqrt(IC_XU6+0.25)+3*U_T')
    i1.defPar('VIN_min', 'V_GS_XU1+Vsat_XU6')
    #Area
    i1.defPar('AreaBcm', 'Wbcm*Lbcm')

    params2html(i1.circuit)
    i1.setDataType('params')
    i1.stepOn()
    i1.setStepVar('Wbcm')
    i1.setStepStart(500e-9)
    i1.setStepStop(30e-6)
    i1.setStepNum(10)
    i1.setStepMethod('log')
    result = i1.execute()
    Sweep1 = plotSweep('SizingBias1', '$CMRR(L_{bcm})$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lbcm', xUnits = 'm', yVar = 'CMRR', yUnits = '', funcType = 'param', show = False)
    fig2html(Sweep1, 600)
    Sweep2 =  plotSweepEdited('SizingBias2', '$Swing: Vin_{min}$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lbcm', xUnits = 'm', yVar = 'VIN_min', yUnits = 'V', funcType = 'param', show = False)
    fig2html(Sweep2, 600)
    text2html('Target: '+str(rangeADC/2))
    Sweep3 = plotSweep('SizingBias3', '$AreaBcm(L_{bcm})$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lbcm', xUnits = 'm', yVar = 'AreaBcm', yUnits = 'm2', funcType = 'param', show = False)
    fig2html(Sweep3, 600)
    i1.stepOff()

##############################################################################
##############################################################################
##############################################################################
##############################################################################
############# OUTPUT STAGE ADJUSTMENT## ######################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################


if Name=="TwoStages" or Name=="FreqComp_Both_NoBias":

        head2html('Trade-offs')
        i1.defPar('Av_stg2ideal', 'g_m_XU4*(1/(g_o_XU4))') #considering circuit 5 (it would make sense to create a circuit with the biasing transistor i guess)
        i1.defPar('Vsat_XU4','2*U_T*sqrt(IC_XU4+0.25)+3*U_T')
        i1.defPar('Vout_max', 'Vdd-Vsat_XU4')
        i1.defPar('AreaO', 'Wout*Lout')
        params2html(i1.circuit)
        i1.setDataType('params')
        i1.stepOn()
        i1.setStepVar('Wout')
        i1.setStepStart(500e-9)
        i1.setStepStop(60e-6)
        i1.setStepNum(10)
        i1.setStepMethod('log')
        result = i1.execute()
        Sweep1 = plotSweep('SizingOut1', 'Gain of the second stage', result, 350e-9, 30e-6, 50, axisType = 'lin', sweepVar= 'Lout', xUnits = 'm', yVar = 'Av_stg2ideal', yUnits = '', funcType = 'param', show = False)
        fig2html(Sweep1, 600)
        Sweep2 = plotSweepEdited('SizingOut2', '$Swing: Vout_{max}$', result, 350e-9, 30e-6, 50, axisType = 'lin', sweepVar= 'Lout', xUnits = 'm', yVar = 'Vout_max', yUnits = 'V', funcType = 'param', show = False)
        fig2html(Sweep2, 600)
        i1.stepOff()




if (Name=="Biasing2" or Name=="FinalCircuit" or Name=="ParasiticCM" or Name=="ParasiticNET108" or Name=="ParasiticOUT1" or Name=="AllParasitics"):
        head2html('Trade-offs')
        i1.defPar('Av_stg2ideal', 'g_m_XU4*(1/(g_o_XU4))')
        i1.defPar('Av_stg2', 'g_m_XU4*(1/(g_o_XU4+g_o_XU5))')
        i1.defPar('GEF_out', 'Av_stg2/Av_stg2ideal')
        i1.defPar('Vsat_XU4','2*U_T*sqrt(IC_XU4+0.25)+3*U_T')
        i1.defPar('Vout_max', 'Vdd-Vsat_XU4')
        i1.defPar('Vsat_XU5','2*U_T*sqrt(IC_XU5+0.25)+3*U_T')
        i1.defPar('Vout_min', 'Vsat_XU5')
        i1.defPar('AreaB', 'Wb*Lb')
        i1.defPar('Vov_XU1', 'V_GS_XU1-Vth_N18')
        i1.defPar('Vov_XU2', 'V_GS_XU2-Vth_P18')
        i1.defPar('Vov_XU3', 'V_GS_XU3-Vth_P18')
        i1.defPar('Vov_XU4', 'V_GS_XU4-Vth_P18')
        i1.defPar('Vov_XU5', 'V_GS_XU5-Vth_N18')
        i1.defPar('Vov_XU6', 'V_GS_XU6-Vth_N18')
        i1.defPar('VeffSI_XU1','IDas/g_m_XU1')
        i1.defPar('VeffSI_XU2','2*IDas/g_m_XU2')
        i1.defPar('VeffSI_XU3','2*IDas/g_m_XU3')
        i1.defPar('VeffSI_XU4','2*IDout/g_m_XU4')
        i1.defPar('VeffSI_XU5','2*IDout/g_m_XU5')
        i1.defPar('VeffSI_XU6','2*2*IDas/g_m_XU6')
        i1.defPar('Vsat_XU1','(2*U_T*sqrt(IC_XU1+0.25)+3*U_T)/2')
        i1.defPar('Vsat_XU2','2*U_T*sqrt(IC_XU2+0.25)+3*U_T')
        i1.defPar('Vsat_XU3','2*U_T*sqrt(IC_XU3+0.25)+3*U_T')
        i1.defPar('Vsat_XU6','2*U_T*sqrt(IC_XU6+0.25)+3*U_T')

        params2html(i1.circuit)
        i1.setDataType('params')
        i1.stepOn()
        i1.setStepVar('Wb')
        i1.setStepStart(500e-9)
        i1.setStepStop(60e-6)
        i1.setStepNum(10)
        i1.setStepMethod('log')
        result = i1.execute()
        Sweep1 = plotSweep('SizingBias1', '$Avstg2(L_{b})$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lb', xUnits = 'm', yVar = 'GEF_out', yUnits = '', funcType = 'param', show = False)
        fig2html(Sweep1, 600)
        Sweep2 = plotSweepEdited('SizingBias2', '$Swing: Vout_{min}$', result, 350e-9, 10e-6, 50, axisType = 'lin', sweepVar= 'Lb', xUnits = 'm', yVar = 'Vout_min', yUnits = 'V', funcType = 'param', show = False)
        fig2html(Sweep2, 600)
        i1.stepOff()

        i1.delPar('Av_stg2ideal')
        i1.delPar('Av_stg2')
        i1.delPar('GEF_out')
        i1.delPar('Vsat_XU4')
        i1.delPar('Vout_max')
        i1.delPar('Vsat_XU5')
        i1.delPar('Vout_min')
        i1.delPar('AreaB')
        htmlPage('Transistor sizing estimation')
        columns = ('W', 'L', 'ID', 'IC', 'gm', 'ro', 'suggested fingers')
        rows = ('U1L','U1R','U2','U3','U4','U5','U6')
        datatable=np.empty((len(rows),len(columns)))
        datatable[0][0]=i1.getParValue('Was')
        datatable[0][1]=i1.getParValue('L')
        datatable[0][2]=i1.getParValue('IDas')
        datatable[0][3]=sp.N(i1.getParValue('IC_XU1'),2)
        datatable[0][4]=sp.N(i1.getParValue('g_m_XU1'),2)
        datatable[0][5]=sp.N(1/i1.getParValue('g_o_XU1'),2)
        datatable[0][6]=sp.N((datatable[0][0]*Rs_N18*datatable[0][4]*(1+datatable[0][3])/(3*datatable[0][1]*(1/2+2/3*datatable[0][3])*i1.getParValue('N_s_N18'))*(GteChnNoiseFactor**2))**(1/2),0)
        #datatable[0][6]=datatable[0][4]/datatable[0][2]

        datatable[1][0]=datatable[0][0]
        datatable[1][1]=datatable[0][1]
        datatable[1][2]=datatable[0][2]
        datatable[1][3]=datatable[0][3]
        datatable[1][4]=datatable[0][4]
        datatable[1][5]=datatable[0][5]
        datatable[1][6]=sp.N((datatable[0][0]*Rs_N18*datatable[0][4]*(1+datatable[0][3])/(3*datatable[0][1]*(1/2+2/3*datatable[0][3])*i1.getParValue('N_s_N18'))*(GteChnNoiseFactor**2))**(1/2),0)
        #datatable[1][6]=datatable[0][6]

        datatable[2][0]=i1.getParValue('Wcm')
        datatable[2][1]=i1.getParValue('Lcm')
        datatable[2][2]=-i1.getParValue('IDasPMOS')
        datatable[2][3]=sp.N(i1.getParValue('IC_XU2'),2)
        datatable[2][4]=sp.N(i1.getParValue('g_m_XU2'),2)
        datatable[2][5]=sp.N(1/i1.getParValue('g_o_XU2'),2)
        datatable[2][6]=sp.N((datatable[2][0]*Rs_P18*datatable[2][4]*(1+datatable[2][3])/(3*datatable[2][1]*(1/2+2/3*datatable[2][3])*i1.getParValue('N_s_P18'))*(GteChnNoiseFactor**2))**(1/2),0)
        #datatable[2][6]=datatable[2][4]/datatable[2][2]

        datatable[3][0]=i1.getParValue('Wcm')
        datatable[3][1]=i1.getParValue('Lcm')
        datatable[3][2]=-i1.getParValue('IDasPMOS')
        datatable[3][3]=sp.N(i1.getParValue('IC_XU3'),2)
        datatable[3][4]=sp.N(i1.getParValue('g_m_XU3'),2)
        datatable[3][5]=sp.N(1/i1.getParValue('g_o_XU3'),2)
        datatable[3][6]=sp.N((datatable[3][0]*Rs_P18*datatable[3][4]*(1+datatable[3][3])/(3*datatable[3][1]*(1/2+2/3*datatable[3][3])*i1.getParValue('N_s_P18'))*(GteChnNoiseFactor**2))**(1/2),0)
        #datatable[3][6]=datatable[3][4]/datatable[3][2]

        datatable[4][0]=i1.getParValue('Wout')
        datatable[4][1]=i1.getParValue('Lout')
        datatable[4][2]=-i1.getParValue('IDoutPMOS')
        datatable[4][3]=sp.N(i1.getParValue('IC_XU4'),2)
        datatable[4][4]=sp.N(i1.getParValue('g_m_XU4'),2)
        datatable[4][5]=sp.N(1/i1.getParValue('g_o_XU4'),2)
        datatable[4][6]=sp.N((datatable[4][0]*Rs_P18*datatable[4][4]*(1+datatable[4][3])/(3*datatable[4][1]*(1/2+2/3*datatable[4][3])*i1.getParValue('N_s_P18'))*(GteChnNoiseFactor**2))**(1/2),1)
        #datatable[4][6]=datatable[4][4]/datatable[4][2]

        datatable[5][0]=i1.getParValue('Wb')
        datatable[5][1]=i1.getParValue('Lb')
        datatable[5][2]=i1.getParValue('IDout')
        datatable[5][3]=sp.N(i1.getParValue('IC_XU5'),2)
        datatable[5][4]=sp.N(i1.getParValue('g_m_XU5'),2)
        datatable[5][5]=sp.N(1/i1.getParValue('g_o_XU5'),2)
        datatable[5][6]=sp.N((datatable[5][0]*Rs_N18*datatable[5][4]*(1+datatable[5][3])/(3*datatable[5][1]*(1/2+2/3*datatable[5][3])*i1.getParValue('N_s_N18'))*(GteChnNoiseFactor**2))**(1/2),1)
        #datatable[5][6]=datatable[5][4]/datatable[5][2]

        datatable[6][0]=i1.getParValue('Wbcm')
        datatable[6][1]=i1.getParValue('Lbcm')
        datatable[6][2]=2*i1.getParValue('IDas')
        datatable[6][3]=sp.N(i1.getParValue('IC_XU6'),2)
        datatable[6][4]=sp.N(i1.getParValue('g_m_XU6'),2)
        datatable[6][5]=sp.N(1/i1.getParValue('g_o_XU6'),2)
        datatable[6][6]=sp.N((datatable[6][0]*Rs_N18*datatable[6][4]*(1+datatable[6][3])/(3*datatable[6][1]*(1/2+2/3*datatable[6][3])*i1.getParValue('N_s_N18'))*(GteChnNoiseFactor**2))**(1/2),1)
        #datatable[6][6]=datatable[6][4]/datatable[6][2]

        fig, ax = plt.subplots()
        ax.set_axis_off()
        table = ax.table(
        cellText = datatable,
        rowLabels = rows,
        colLabels = columns,
        rowColours =["mediumturquoise"] * (len(rows)+1), #
        colColours =["mediumturquoise"] * (len(columns)+1), #
        cellLoc ='center',
        loc ='upper left')
        ax.set_title('Transistor sizing',
                 fontweight ="bold")
        fig.savefig('./img/TransistorSizingTable.svg',dpi=600)
        img2html('TransistorSizingTable.svg',600)

        htmlPage('OpAmp performance summary')
        #Gain
        i1.defPar('Av_stg2', 'g_m_XU4*(1/(g_o_XU4+g_o_XU5))')
        i1.defPar('Av_stg1', '(g_m_XU1)*(1/((g_o_XU1)+g_o_XU3))*(2*g_m_XU3*(1/g_o_XU3)+1)/(2*(g_m_XU3*(1/g_o_XU3)+1))')
        Av1=i1.getParValue('Av_stg1')
        Av2=i1.getParValue('Av_stg2')
        #CMIR
        i1.defPar('Vsat_XU1','(2*U_T*sqrt(IC_XU1+0.25)+3*U_T)/2')
        i1.defPar('MaxIN', 'Vdd+V_GS_XU3-Vsat_XU1+V_GS_XU1')
        i1.defPar('Vsat_XU6','2*U_T*sqrt(IC_XU6+0.25)+3*U_T')
        i1.defPar('VIN_min', 'V_GS_XU1+Vsat_XU6')
        CMIRH=sp.N(i1.getParValue('MaxIN'),3)
        CMIRL=sp.N(i1.getParValue('VIN_min'),3)
        text2html('Vin range: from '+str(CMIRL)+' to '+str(CMIRH)+' V')
        #Output bias range
        i1.defPar('Vsat_XU4','2*U_T*sqrt(IC_XU4+0.25)+3*U_T')
        i1.defPar('Vout_max', 'Vdd-Vsat_XU4')
        i1.defPar('Vsat_XU5','2*U_T*sqrt(IC_XU5+0.25)+3*U_T')
        i1.defPar('Vout_min', 'Vsat_XU5')
        VBOUTH=sp.N(i1.getParValue('Vout_max'),3)
        VBOUTL=sp.N(i1.getParValue('Vout_min'),3)
        text2html('Vout range: from '+str(VBOUTH)+' to '+str(VBOUTL)+' V')
        #PSRR
        i1.defPar('PSRR','g_m_XU1/(g_o_XU1+g_o_XU3)')
        PSRR=i1.getParValue('PSRR')
        text2html('PSRR: '+str(sp.N(20*log(PSRR,10),3))+' dB')
        i1.defPar('PSRR2','g_m_XU4/g_o_XU4')
        PSRR2=i1.getParValue('PSRR2')
        text2html('PSRR2: '+str(sp.N(20*log(PSRR2,10),3))+' dB, effective: '+str(sp.N(20*log(PSRR2*Av1,10),3)))
        #CMRR
        i1.defPar('CMRR', '(1+2*(g_m_XU1)*(1/g_o_XU6))*g_m_XU3*(1/(g_o_XU1+g_o_XU3))')
        CMRR=i1.getParValue('CMRR')
        text2html('CMRR: '+str(sp.N(20*log(CMRR,10),3))+' dB')
        #Power
        Power_tot=sp.N((datatable[0][2]*2+datatable[5][2])*Vdd,3)
        text2html('Power: '+str(Power_tot)+' W')
        #Area
        AreaDP=sp.N(datatable[0][0]*datatable[0][1]*2,3)
        AreaCM=sp.N(datatable[2][0]*datatable[2][1]*2,3)
        AreaCS=sp.N(datatable[4][0]*datatable[4][1],3)
        AreaB2=sp.N(datatable[5][0]*datatable[5][1],3)
        AreaB1=sp.N(datatable[6][0]*datatable[6][1],3)
        Area_tot=sp.N(AreaDP+AreaCM+AreaCS+AreaB2+AreaB1+CompensationArea,3)
        text2html('Area: '+str(Area_tot)+' m2')
        text2html('     Diff.Pair: '+str(AreaDP)+' m2')
        text2html('     Curr.mirr: '+str(AreaCM)+' m2')
        text2html('     Cs.output: '+str(AreaCS)+' m2')
        text2html('     Bias.outp: '+str(AreaB2)+' m2')
        text2html('     Bias.inpu: '+str(AreaB1)+' m2')
        text2html('     Compensation: '+str(CompensationArea)+' m2')
        #OpAmp noise (traditional)
        i1.defPar('NEF', '1+(g_m_XU3/(g_m_XU1))') #"Noise excess factor"
        NEF=sp.N(i1.getParValue('NEF'),3)
        text2html('NEF: '+str(NEF)+' increment in respect to the noise in design step 2')
        #Estimated GAIN
        dcL=float(sp.N(loopgainpz.DCvalue, ini.disp))
        text2html('TIA Open loop gain: '+str(sp.N(20*log(-dcL,10),3))+' dB')
        text2html('OpAmp gain: '+str(sp.N(20*log(Av1*Av2,10),3))+' dB')
        text2html(' 1st stage: '+str(sp.N(20*log(Av1,10),3))+' dB')
        text2html(' 2nd stage: '+str(sp.N(20*log(Av2,10),3))+' dB')
        text2html('Rps='+str(format(i1.getParValue('Rps'),'.2E'))+'$\Omega$')
        text2html('Cps='+str(format(i1.getParValue('Cps'),'.2E'))+'$\Omega$')
        text2html('Rph1='+str(format(i1.getParValue('Rph1'),'.2E'))+'$\Omega$')

        #Estimated bandwidth (worst case, max_gain)
        i1.defPar('gain', max_gain)
        i1.setSimType('numeric') #numeric simulation to substitute all variables with their values
        i1.setSource('I1')       #the source is the membrane current source
        i1.setDetector('V_out')  #the detector is the output (input of the ADC)
        i1.setLGref(reference) #reference source in the asymptotic model
        i1.setGainType('loopgain')
        i1.setDataType('laplace')
        loopgain1 = i1.execute()
        servoData = findServoBandwidth(loopgain1.laplace)      #calculate the frequency cut off and the midband value of the loopgain among others
        f_h = servoData['lpf']
        text2html('OpAmp BW (worst case max_gain): '+str(sp.N(f_h,3))+' Hz')

        i1.defPar('gain', m_gain)
        loopgain2 = i1.execute()
        servoData = findServoBandwidth(loopgain2.laplace)      #calculate the frequency cut off and the midband value of the loopgain among others
        f_h = servoData['lpf']
        text2html('OpAmp BW (m_gain): '+str(sp.N(f_h,3))+' Hz')

        i1.defPar('gain', min_gain)
        loopgain3 = i1.execute()
        servoData = findServoBandwidth(loopgain3.laplace)      #calculate the frequency cut off and the midband value of the loopgain among others
        f_h = servoData['lpf']
        text2html('OpAmp BW (min_gain): '+str(sp.N(f_h,3))+' Hz')

        #PLOT BODE
        ini.plotFontSize=18
        head2html('Max gain transfers')
        i1.defPar('gain', max_gain) #to check the direct transfer issue
        i1.setGainType('gain')
        i1.setDataType('laplace')
        gain = i1.execute()
        i1.setGainType('asymptotic')
        i1.setDataType('laplace')
        asymptotic = i1.execute()
        i1.setGainType('servo')
        i1.setDataType('laplace')
        servo = i1.execute()
        i1.setGainType('direct')
        i1.setDataType('laplace')
        direct = i1.execute()
        head2html('Bode plots')
        result = [asymptotic, gain, loopgain1, servo, direct]
        figdBmag = plotSweep('dBmagFinal', 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
        figPhase = plotSweep('phaseFinal','Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
        fig2html(figdBmag, 800)
        fig2html(figPhase, 800)
        print("Bode plots obtained")
        #params2html(i1.circuit)

        #PLOT BODE (Checking m_gain case)
        ini.plotFontSize=18
        head2html('M gain transfers')
        i1.defPar('gain', m_gain) #to check the direct transfer issue
        i1.setGainType('gain')
        i1.setDataType('laplace')
        gain = i1.execute()
        i1.setGainType('asymptotic')
        i1.setDataType('laplace')
        asymptotic = i1.execute()
        i1.setGainType('servo')
        i1.setDataType('laplace')
        servo = i1.execute()
        i1.setGainType('direct')
        i1.setDataType('laplace')
        direct = i1.execute()
        head2html('Bode plots')
        result = [asymptotic, gain, loopgain2, servo, direct]
        figdBmag2 = plotSweep('dBmagFinal2', 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
        figPhase2 = plotSweep('phaseFinal2','Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
        fig2html(figdBmag2, 800)
        fig2html(figPhase2, 800)
        print("Bode plots obtained")

        i1.setGainType('gain')
        i1.setDataType('pz')
        gainpz=i1.execute()
        pz2html(gainpz)

        #PLOT BODE (Checking m_gain case)
        ini.plotFontSize=18
        head2html('Min gain transfers')
        i1.defPar('gain', min_gain) #to check the direct transfer issue
        i1.setGainType('gain')
        i1.setDataType('laplace')
        gain = i1.execute()
        i1.setGainType('asymptotic')
        i1.setDataType('laplace')
        asymptotic = i1.execute()
        i1.setGainType('servo')
        i1.setDataType('laplace')
        servo = i1.execute()
        i1.setGainType('direct')
        i1.setDataType('laplace')
        direct = i1.execute()
        head2html('Bode plots')
        result = [asymptotic, gain, loopgain3, servo, direct]
        figdBmag3 = plotSweep('dBmagFinal3', 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
        figPhase3 = plotSweep('phaseFinal3','Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
        fig2html(figdBmag3, 800)
        fig2html(figPhase3, 800)
        print("Bode plots obtained")
