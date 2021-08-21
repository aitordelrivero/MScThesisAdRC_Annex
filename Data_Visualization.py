##!/usr/bin/env python2
# -*- coding: utf-8 -*-

from SLiCAP import *
from sympy.plotting import plot3d
from scipy.optimize import fminbound
from math import ceil, prod, log
from itertools import product
from scipy.integrate import quad

def AssignProcessParameters(self):
    """
    Sets simplified EKV parameters

    """
    #######???????????????????????????????
    #COPY HERE PROCESS PARAMETERS AVAILABLE IN THE WIKI
    #ELSE, THE SCRIPT USES DEFAULT CMOS18 EKV PARAMETERS FROM BINKLEY
    #######???????????????????????????????
    pass

def Cadencecsv2traces(csvFile, absx = False, logx = False, absy = False, logy = False, limiter=None):
    """
    Note: is similar to "csv2traces" in https://github.com/Lenty/SLiCAP_python.
    The csv delimiter is modifed to understand the multi-sweep Cadence csv generation terminology.
    There are additional functions to process the csv data (abs, log)
    "limiter" fixes the bug in transient simulations, where the last pairs (X,Y) are not well represented

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
    for label in labels:
            traceDict[label].xData = np.array(traceDict[label].xData)
            traceDict[label].yData = np.array(traceDict[label].yData)
            traceDict[label].label = label
    return traceDict

def LoadFinalCircuit(Parasitics=False,LoopgainMag=True,LoopgainPhs=True):
    """
    Returns the Loop Gain symbolic analysis estimation of the final circuit

    Select Parasitcs=True for symbolic estimations considering the parasitics after extraction in Cadence.

    """
    #Frequency limits to be plot
    f_check_max=1e10
    f_check_min=1e-3

    # RESULTS COMING FROM OTHER SCRIPTS
    Wopt=30e-6                      #Optimum width for the input stage (from Script "Noise study")
    Iopt=2e-6                       #Optimum current for the input stage (from Script "Noise study")
    L = 2e-6                        #Optimum length for the input stage (from Script "Noise study")
    Wout=5e-6                       #Output stage transistor width (From Script "Stages Study")
    Lout=350e-9                     #Output stage transistor length (From Script "Stages Study")
    IDout=5e-6                      #Output stage transistor drain current bias (From Script "Stages Study")
    Lcm=5e-6                        #1st stage Top biasing transistor length (From Script "Stages Study")
    Wcm=7.5e-6                      #1st stage Top biasing transistor width (From Script "Stages Study")
    Wb=5e-6                         #2nd stage biasing transistor width (From Script "Stages Study")
    Lb=1.5e-6                       #2nd stage biasing transistor length (From Script "Stages Study")
    Wbcm=8e-6                       #1st stage Bottom biasing transistor width (From Script "Stages Study")
    Lbcm=1.5e-6                     #1st stage  Bottom biasing transistor length (From Script "Stages Study")

    CpsDefault='5p'                 #Pole splitting capacitance  (from Script Freq. Compensation Study)
    RpsDefault=200000               #Pole splitting resistance (from Script Freq. Compensation Study)
    Rph1Default=4000                #Phantom zero splitting resistance (from Script Freq. Compensation Study)

    CompensationArea=26.5e-6**2+(4+1)*350e-9*10e-6+1e-6*10e-6 #Cps+RpsPMOS+RpsNMOS+Rph1 areas

    CparX=10.5e-15                                            #Parasitic at node X (from Cadence)
    CparC=25.03e-15                                           #Parasitic at node C (from Cadence)
    CparY=29.73e-15                                           #Parasitic at node Y (from Cadence)


    IG = 0                          #Gate current is neglectible
    Ce = 5e-12                      #Electrode impedance
    Ree = 1500e9                    #Electrode resistance
    Rs = 100e6                      #Solution resistance
    f_min=1e-4                      #Minimum interest frequency (aprox. DC)
    f_max = 10000                   #Maximum interest frequency
    max_gain=20e9                   #Maxiumum gain of the TIA
    min_gain=2e6                    #Minimum gain of the TIA
    m_gain=200e6                    #Intermediate gain of the TIA

    Cadc=5e-12                      #Equivalent input capacitance of the ADC
    Cgm=0                           #Equivalent capacitance of the feedback network
    Cload=Cadc+Cgm
    res=10                          #Resolution (x-bit) of the ADC
    vRange=2                        #Voltage range of the ADC
    f_sampling=20e3                 #Sampling frequency of the ADC

    Vdd=3.3                         #Power supply


    ngtbc=3 #Number of gains to be checked


    #Selection of estimations with or without parasitics
    if not Parasitics:
        circuit="FinalCircuit"
    if Parasitics:
        circuit="AllParasitics"
    Name=circuit
    test=False

    makeNetlist(Name+'.asc',Name)       #create netlist
    i1=instruction()                    #create instruction object
    i1.setCircuit(Name+'.cir')          #load circuit
    print("-------------------")
    print("CIRCUIT SET: "+Name)
    AssignProcessParameters(i1)

    if (Name=="CIRCUIT1" or Name=="CIRCUIT2" or Name=="CIRCUIT3" or Name=="CIRCUIT3b" or Name=="CIRCUIT5" or Name=="CIRCUIT6" or Name=="FinalCircuit" or Name=="ParasiticCM" or Name=="ParasiticNET108" or Name=="ParasiticOUT1" or Name=="AllParasitics"):
        reference='Gm_M1_XU1' #select loop reference source of the asymptotic gain model
        i1.defPar('c_dg_XU1', 0) #eliminate local feedback loop around the reference
        print('loop gain reference assigned, and local loop in main signal path deleted')
    if (Name=="CIRCUIT3" or Name=="CIRCUIT3b" or Name=="CIRCUIT5" or Name=="CIRCUIT6" or Name=="FinalCircuit" or Name=="ParasiticCM" or Name=="ParasiticNET108" or Name=="ParasiticOUT1" or Name=="AllParasitics"):
        i1.defPar('c_dg_XU3', 0) #eliminate local feedback loop around the reference
        i1.defPar('c_dg_XU2', 0) #eliminate local feedback loop around the reference
        print('local loop in secondary signal path deleted')
    else:
        print('loop gain reference was not assigned')
    try: #Assigns variables to the input stage elements
        i1.defPar('W', Wopt)
        i1.defPar('ID', Iopt)
        i1.defPar('L', L)
        print('single-stage parameters assigned')
    except:
        print('single-stage parameters not needed')
    try: #Assigns variables to the output stage elements
        i1.defPar('IDoutPMOS', -IDout)
        i1.defPar('Wout', Wout)
        i1.defPar('Lout', Lout)
        i1.defPar('Vdd', Vdd)
        i1.defPar('range_ADC', range)
        print('output parameters assigned')
    except:
        print('two-stage ouput parameters not needed')
    try: #Assigns variables to the balanced input stage elements
        i1.defPar('IDas', Iopt*2)
        i1.defPar('Was', Wopt*2)
        i1.defPar('Lin', L)
        print('input balanced CS parameters assigned')
    except:
        print('input antiseries not needed')
    try: #Assigns variables to the input stage top biasing elements
        i1.defPar('Wcm', Wcm)
        i1.defPar('Lcm', Lcm)
        i1.defPar('IDasPMOS', -Iopt*2)
        print('current mirror parameters assigned')
    except:
        print('no current mirrors')
    try: #Assigns variables to the output stage biasing elements
        i1.defPar('Wb', Wb)
        i1.defPar('Lb', Lb)
        i1.defPar('IDout', IDout)
        print('output bias equivalent parameters assigned')
    except:
        print('no output bias equivalent')
    try: #Assigns variables to the input stage bottom biasing elements
        i1.defPar('Wbcm', Wbcm)
        i1.defPar('Lbcm', Lbcm)
        i1.defPar('IDas', Iopt*2)
        print('output bias equivalent parameters assigned')
    except:
        print('no input bias equivalent')
    try:
        i1.defPar('CparX', CparX) #Assigns parasitic extraction values
        print('parasitic CparX set')
    except:
        print('parasitic CparX not modeled')
    try:
        i1.defPar('CparC', CparC) #Assigns parasitic extraction values
        print('parasitic CparC set')
    except:
        print('parasitic CparC not modeled')
    try:
        i1.defPar('CparY', CparY) #Assigns parasitic extraction values
        print('parasitic CparY set')
    except:
        print('parasitic CparY not modeled')


    # Assigns model specefications
    i1.defPar('IG', IG)
    i1.defPar('Ce', Ce)
    i1.defPar('Ree', Ree)
    i1.defPar('Rs', Rs)
    i1.defPar('Cload', Cload)
    print("Parameters assigned")

    # Assigns compensation elements' values
    i1.defPar('Rps', RpsDefault)
    i1.defPar('Cps', CpsDefault)
    i1.defPar('Rph1', Rph1Default)
    i1.defPar('gain', m_gain)

    i1.setSimType('numeric') #numeric simulation to substitute all variables with their values
    i1.setSource('I1')       #the source is the membrane current source
    i1.setDetector('V_out')  #the detector is the output (input of the ADC)
    i1.setLGref(reference) #reference source in the asymptotic model
    i1.setStepVar('gain')
    #Estimates the loop gain of the amplifier
    i1.setGainType('loopgain')
    i1.setDataType('laplace')
    loopgain = i1.execute()
    head2html('Bode plot. Loopgain.')
    result = [loopgain]
    if LoopgainMag:
        SLiCAPPlotSweepMag = plotSweep('LaTEX_ResultSLICAPLoopgainMag_Chapter'+str(CHAPTER), 'Loop gain magnitude (Gain=200MV/A)', result, f_check_min, f_check_max, 100, funcType = 'dBmag', show=False)
    else:
        SLiCAPPlotSweepMag="Void"
    if LoopgainPhs:
        SLiCAPPlotSweepPhase = plotSweep('LaTEX_ResultSLICAPLoopgainPhs_Chapter'+str(CHAPTER), 'Loop gain phase (Gain=200MV/A)', result, f_check_min, f_check_max, 100, funcType = 'phase', show=False)
    else:
        SLiCAPPlotSweepPhase="Void"
    print("Loopgain bode obtained")
    #returns the loop gain
    return f_check_max,f_check_min,SLiCAPPlotSweepMag, SLiCAPPlotSweepPhase

def LoadBasicCircuit(Name, Ce, Ree,Rs):
    i1=instruction()
    i1.setCircuit(Name+'.cir')
    print("Circuit set")
    AssignProcessParameters(i1)
    i1.defPar('IG', 0)
    i1.defPar('Ce', Ce)
    i1.defPar('Ree', Ree)
    i1.defPar('Rs', Rs)
    return i1

def LoadNoiseSymbolic():
    Wopt=30e-6                      #Optimum width for the input stage (from Script "Noise study")
    Iopt=2e-6                       #Optimum current for the input stage (from Script "Noise study")
    L = 2e-6                        #Optimum length for the input stage (from Script "Noise study")

    NEF=1.25                        #Noise Excess Factor, (from Script "Stage study")

    IG = 0                          #Gate current is neglectible
    Ce = 5e-12                      #Electrode impedance
    Ree = 1500e9                    #Electrode resistance
    Rs = 100e6                      #Solution resistance
    f_min=1e-4                      #Minimum interest frequency (aprox. DC)
    f_max = 10000                   #Maximum interest frequency
    max_gain=20e9                   #Maxiumum gain of the TIA
    min_gain=2e6                    #Minimum gain of the TIA
    m_gain=200e6                    #Intermediate gain of the TIA

    res=10                          #Resolution (x-bit) of the ADC
    vRange=2                        #Voltage range of the ADC
    f_sampling=20e3                 #Sampling frequency of the ADC
    noise_budget_vs_feedback=50/100 #Noise budget of the amplifier: splitted 50-50 between the feedback network and the controller.

    #Execution parameters:
    ngtbc=100 #Number of gains to be checked

    #Establish noise target
    LSB=vRange/(2**res)                                             #LSB of the ADC
    target=LSB/(sp.sqrt(12))*(sp.sqrt(noise_budget_vs_feedback))    #quantization noise * budget ratio by the controller

    #Obtain estimated noise
    Name='NoiseController'
    htmlPage('Controller noise optimization')
    i1=LoadBasicCircuit(Name,Ce, Ree,Rs) #loads noise model
    i1.defPar('L', L)                    #assigns paramters
    i1.defPar('W', Wopt)
    i1.defPar('ID', Iopt)
    gain_sel  = sp.Symbol('gain')        #Parameter to sweep
    #Setting simulation config
    i1.setSimType('numeric')   #To avoid unnecesary computation effort by replacing
                               #with numeric values all variables whose value we know.
                               #Applicable to all variables but W, ID and gain.
    i1.setGainType('vi')       #to obtain transfer functions
    i1.setDataType('noise')    #to calculate noise equations from noise sources at the .cir
    i1.setSource('I1_XU1')     #the only noise source in the model is the controller's noise equivalent
    i1.setDetector('V_out')    #the detector where we measure the noise (input of the ADC = output of the amplifier)
    resultA = i1.execute()

    #plot the rms onoise vs gain
    gain_sel_lists=np.geomspace(min_gain, max_gain, num=ngtbc, endpoint=True).tolist()
    X=np.geomspace(min_gain, max_gain, num=ngtbc, endpoint=True).tolist()
    MyYo  = np.zeros([ngtbc])
    # Create a function for numeric substitution of grid point values
    onoise = sp.lambdify([gain_sel], resultA.onoise)
    for i in range(ngtbc):
            onoiseF = sp.lambdify(ini.frequency, onoise(gain_sel_lists[i]))
            # Calculate RMS value for current grid point
            MyYo[i] = np.sqrt(quad(onoiseF, f_min, f_max)[0])*NEF
    ini.defaultColors=['k','#99d28c','#f1be3e','#007188','k','#eb7245','#0066a2','#00a390','#82d7c6']
    pairs=[X, MyYo]
    SLiCAPconverted = {"Symbolic Analysis": pairs}
    GraphNoiseRMS = plot('onoise_RMS_200', 'Output referred noise (RMS DC-10kHz)', 'log', SLiCAPconverted, xName = 'Gain', xUnits = 'V/A', yName = 'onoise', yUnits = 'Vrms',xLim = [] , yLim = [], show = False, )
    ini.defaultColors      = ['r','b','g','c','m','y','k']

    #plot the spectrum for 200MV/A
    head2html('Spectrum plots for gain=200MV/A')
    i1.defPar('gain', 200e6)
    i1.setSource('I1_XU1')     #the only noise source in the model is the controller's noise equivalent
    result = i1.execute()
    i1.setSource('I1')
    result2 = i1.execute()
    result.onoise=result.onoise*NEF
    result2.onoise=result.onoise*NEF
    NoiseSpectrum=plotSweep('onoise_spect_200','Output referred noise sprectrum. Gain=200M.', result, 8e-3, 2e4, 200, funcType='onoise', show=False) #plot the spectrum

    #return the estimated rms onoise and spectrum with SLICAP
    return NoiseSpectrum, GraphNoiseRMS

t1 = time()

prj = initProject('Visualization')
print("Project inititated")
ini.MaximaTimeOut = 600         #Some extra time to allow the computer to calculate integrals.

#GRAPH SELECTION PARAMETERS (EXECUTION PARAMETERS)
CHAPTER=51                                      #For selecting the loop gains. 42 Controller with ideal resistor, 44 Post-layout controller with ideal resistor, 51 Prelayout full, 52 Postlayout full, 'All' for all together
LoopgainMag=True                                #PLot loop gain magnitude comparison between Cadence simulations and SLiCAP estimations
LoopgainPhs=True                                #PLot loop gain phase comparison between Cadence simulations and SLiCAP estimations
NoiseSpectrum=False                             #PLot noise spectrum comparison between Cadence simulations and SLiCAP estimations
NoiseRMSResult=False                            #PLot onoise rms vs gain comparison between Cadence simulations and SLiCAP estimations
P_Res_Ideal_4vs8_ResponsesDCAC=False            #PLot AC and DC responses of the pseudo-resistor using 4 and 8 cells comparison between Cadence simulations and SLiCAP estimations
PseudoResistor_AC_RealvsIdeal=False             #PLot AC response ideal pseudo resistor and with biasing comparison between Cadence simulations and SLiCAP estimations
PseudoResistor_DC_RealvsIdeal=False             #PLot DC response ideal pseudo resistor and with biasing comparison between Cadence simulations and SLiCAP estimations
ACLevels=False                                  #PLot AC response Pre and Post layout for different gain levels comparison between Cadence simulations and SLiCAP estimations
PseudoresistorBiasing=False                     #PLot biasing Vc generation versus Ic for High-resistance and Low-resistance modes.
ResistanceSelection=False                       #PLot loop gain phase of the amplifier using a a resistor or a parallel of triode transistors to implement Rps
Transient=False                                 #PLot transient (voltage clamp + membrane current stimuli + output voltage)
NoiseChapter51pre=False                         #PLot noise spectrum and rms onoise vs gain for the whole system (pre-layout) and comparison with estimated controller contribution
NoiseChapter51post=False                        #PLot noise spectrum and rms onoise vs gain for the whole system (post-layout) and comparison with estimated controller contribution
gainAC=False                                    #PLot transimpedance gain from membrane to load vs frequency
servoCompare=False                              #Plot pre and post layout servo functions
PSRRiPlot=False                                 #PLot PSRR referred to the membrane current
PSRRvPlot=False                                 #PLot PSRR referred to the clamping voltage
Driving=False                                   #PLot the results of the testbench to determine the ouput capabilities

if (LoopgainMag or LoopgainPhs):
    f_check_max,f_check_min,SLiCAPPlotSweepMag, SLiCAPPlotSweepPhase = LoadFinalCircuit(False,LoopgainMag,LoopgainPhs)

if (NoiseSpectrum or NoiseRMSResult):
    SLiCAPPlotNoise,SLICAPNoiseRMS=LoadNoiseSymbolic()
if NoiseSpectrum:
    ini.plotFontSize=18
    htmlPage("Noise")
    Trace1 = csv2traces('onoise_Aitor.csv')
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
         if key in ['onoise (WRps=1e-06) Y']:
             Trace1sel[key]=Trace1[key]
             Trace1sel[key].label=Trace1sel[key].label.replace("onoise (WRps=1e-06) Y","Controller - Gain 200M")
             print(str(Trace1sel[key].label))
    TracePlot={}
    keys=list(Trace1sel.keys())
    selectionList=['loopgain','Controller - Gain 200M']
    tag1='Symbolic analysis'
    color1='k'
    tag2='Cadence'
    color2='#00a390'
    ini.defaultColors      = [color1,color2]
    for key in keys:
        if Trace1sel[key].label in selectionList:
            TracePlot[key]=Trace1sel[key]
            TracePlot[key].label=TracePlot[key].label.replace("Controller - Gain 200M",tag2)
    SLiCAPPlotNoise.axes[0][0].traces[0].label=SLiCAPPlotNoise.axes[0][0].traces[0].label.replace("onoise",tag1)
    traces2fig(TracePlot,SLiCAPPlotNoise)
    SLiCAPPlotNoise.plot()
    fig2html(SLiCAPPlotNoise, 600)
    del Trace1, Trace1sel, TracePlot
    ini.defaultColors      = ['r','b','g','c','m','y','k']


if NoiseRMSResult:
    ini.plotFontSize=18
    color1='k'
    color2='#00a390'
    ini.defaultColors      = [color1,color2]
    Trace1 = csv2traces('RMSNoiseController.csv')
    Trace1sel={}
    keys=list(Trace1.keys())
    print(keys)
    for key in keys:
        if key in ['Vrms_onoise (Vclamp=1.65) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Vrms_onoise (Vclamp=1.65) Y","Cadence")
            print(str(Trace1sel[key].label))
    traces2fig(Trace1sel,SLICAPNoiseRMS)
    SLICAPNoiseRMS.plot()
    fig2html(SLICAPNoiseRMS, 600)
    ini.defaultColors      = ['r','b','g','c','m','y','k']

if LoopgainMag:
    ini.plotFontSize=18
    Trace1 = Cadencecsv2traces('LoopGainControllerMag.csv')
    Trace2 = Cadencecsv2traces('LoopgainMagPrelayout3.csv')
    Trace3 = Cadencecsv2traces('LoopGainMagPrePostLayoutComparisonOpAmp.csv')
    Trace4 = Cadencecsv2traces('LoopGainMagPrePostLayoutComparisonWholeCircuit.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,Rfb=2e+08) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain dB20 (Vclamp=1.65,Rfb=2e+08) Y","Prelayout. Ideal resistor.")
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("Loop Gain dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y","Prelayout. Pseudoresistor.")
    Trace3sel={}
    keys=list(Trace3.keys())
    for key in keys:
        if key in ['Loop Gain db20 Postlayout (Vclamp=1.65,Rfb=2e+08) Y','Loop Gain db20 Prelayout (Vclamp=1.65,Rfb=2e+08) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain db20 Postlayout (Vclamp=1.65,Rfb=2e+08) Y","Postlayout. Ideal resistor.")
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain db20 Prelayout (Vclamp=1.65,Rfb=2e+08) Y","Prelayout. Ideal resistor. Blackbox.")
    Trace4sel={}
    keys=list(Trace4.keys())
    for key in keys:
        if key in ['Loop Gain dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y','Loop Gain dB20 Postlayout (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain dB20 Postlayout (Vclamp=1.65,A=0,Idac=9e-08) Y","Postlayout. Whole circuit.")
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y","Prelayout. Pseudoresistor. Blackbox.")
    Trace1sel.update(Trace2sel)
    Trace1sel.update(Trace3sel)
    Trace1sel.update(Trace4sel)
    print(list(Trace1sel.keys()))

    TracePlot={}
    keys=list(Trace1sel.keys())
    selectionList=['loopgain']
    tag1='Symbolic analysis'
    color1='k'
    tag2='Prelayout. Ideal resistor.'
    color2='#00a390'
    tag3='Postlayout. Ideal resistor.'
    color3='#00b7d3'
    tag4='Prelayout. Pseudoresistor.'
    color4='#eb7245'
    tag5='Postlayout. Whole circuit.'
    color5='#f1be3e'
    if CHAPTER==42:
        selectionList=['loopgain','Prelayout. Ideal resistor.']
        tag1="Symbolic analysis"
        tag2="Cadence simulation"
        ini.defaultColors      = [color1,color2]
    if CHAPTER==44:
        selectionList=['loopgain','Prelayout. Ideal resistor.','Postlayout. Ideal resistor.']
        tag1="Symbolic analysis"
        tag2="Cadence - Prelayout"
        tag3="Cadence - Postlayout"
        ini.defaultColors      = [color1,color2,color3]
    if CHAPTER==51:
        selectionList=['loopgain','Prelayout. Ideal resistor.','Prelayout. Pseudoresistor.']
        tag1="Symbolic analysis"
        tag2="Cadence w/o pseudoresistor"
        tag4="Cadence w pseudoresistor"
        ini.defaultColors      = [color1,color2,color4]
    if CHAPTER==52:
        selectionList=['loopgain','Postlayout. Ideal resistor.','Postlayout. Whole circuit.']
        tag1="Symbolic analysis"
        tag3="Cadence w/o pseudoresistor"
        tag5="Cadence w pseudoresistor"
        ini.defaultColors      = [color1,color3,color5]
    if (CHAPTER==42 or CHAPTER==44 or CHAPTER==51 or CHAPTER==52):
        for key in keys:
            if Trace1sel[key].label in selectionList:
                TracePlot[key]=Trace1sel[key]
                TracePlot[key].label=TracePlot[key].label.replace("Prelayout. Ideal resistor.",tag2)
                TracePlot[key].label=TracePlot[key].label.replace("Postlayout. Ideal resistor.",tag3)
                TracePlot[key].label=TracePlot[key].label.replace("Prelayout. Pseudoresistor.",tag4)
                TracePlot[key].label=TracePlot[key].label.replace("Postlayout. Whole circuit.",tag5)
    else:
        TracePlot=Trace1sel
        ini.defaultColors      = [color1,color2,color4,color3,'gray10','gray10',color5]
    SLiCAPPlotSweepMag.axes[0][0].traces[0].label=SLiCAPPlotSweepMag.axes[0][0].traces[0].label.replace("loopgain",tag1)
    traces2fig(TracePlot,SLiCAPPlotSweepMag)
    SLiCAPPlotSweepMag.plot()
    fig2html(SLiCAPPlotSweepMag, 600)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel, TracePlot
    ini.defaultColors      = ['r','b','g','c','m','y','k']


if LoopgainPhs:
    ini.plotFontSize=18
    Trace1 = Cadencecsv2traces('LoopGainControllerPhase.csv')
    Trace2 = Cadencecsv2traces('LoopgainPhasePrelayout3.csv')
    Trace3 = Cadencecsv2traces('LoopGainPhasePrePostLayoutComparisonOpAmp.csv')
    Trace4 = Cadencecsv2traces('LoopGainPhasePrePostLayoutComparisonWholeCircuit.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,Rfb=2e+08) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Loop Gain Phase (Vclamp=1.65,Rfb=2e+08) Y","Prelayout. Ideal resistor.")
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("Loop Gain Phase (Vclamp=1.65,A=0,Idac=9e-08) Y","Prelayout. Pseudoresistor.")
    Trace3sel={}
    keys=list(Trace3.keys())
    for key in keys:
        if key in ['Loop Gain Phase Postlayout (Vclamp=1.65,Rfb=2e+08) Y','Loop Gain Phase Prelayout (Vclamp=1.65,Rfb=2e+08) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain Phase Postlayout (Vclamp=1.65,Rfb=2e+08) Y","Postlayout. Ideal resistor.")
            Trace3sel[key].label=Trace3sel[key].label.replace("Loop Gain Phase Prelayout (Vclamp=1.65,Rfb=2e+08) Y","Prelayout. Ideal resistor. Blackbox.")
    Trace4sel={}
    keys=list(Trace4.keys())
    for key in keys:
        if key in ['Loop Gain Phase (Vclamp=1.65,A=0,Idac=9e-08) Y','Loop Gain Phase Postlayout (Vclamp=1.65,A=0,Idac=9e-08) Y']:
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain Phase Postlayout (Vclamp=1.65,A=0,Idac=9e-08) Y","Postlayout. Whole circuit.")
            Trace4sel[key].label=Trace4sel[key].label.replace("Loop Gain Phase (Vclamp=1.65,A=0,Idac=9e-08) Y","Prelayout. Pseudoresistor. Blackbox.")
    Trace1sel.update(Trace2sel)
    Trace1sel.update(Trace3sel)
    Trace1sel.update(Trace4sel)
    print(list(Trace1sel.keys()))

    TracePlot={}
    keys=list(Trace1sel.keys())
    selectionList=['loopgain']
    tag1='Symbolic analysis'
    color1='k'
    tag2='Prelayout. Ideal resistor.'
    color2='#00a390'
    tag3='Postlayout. Ideal resistor.'
    color3='#00b7d3'
    tag4='Prelayout. Pseudoresistor.'
    color4='#eb7245'
    tag5='Postlayout. Whole circuit.'
    color5='#f1be3e'
    if CHAPTER==42:
        selectionList=['loopgain','Prelayout. Ideal resistor.']
        tag1="Symbolic analysis"
        tag2="Cadence simulation"
        ini.defaultColors      = [color1,color2]
    if CHAPTER==44:
        selectionList=['loopgain','Prelayout. Ideal resistor.','Postlayout. Ideal resistor.']
        tag1="Symbolic analysis"
        tag2="Cadence - Prelayout"
        tag3="Cadence - Postlayout"
        ini.defaultColors      = [color1,color2,color3]
    if CHAPTER==51:
        selectionList=['loopgain','Prelayout. Ideal resistor.','Prelayout. Pseudoresistor.']
        tag1="Symbolic analysis"
        tag2="Cadence w/o pseudoresistor"
        tag4="Cadence w pseudoresistor"
        ini.defaultColors      = [color1,color2,color4]
    if CHAPTER==52:
        selectionList=['loopgain','Postlayout. Ideal resistor.','Postlayout. Whole circuit.']
        tag1="Symbolic analysis"
        tag3="Cadence w/o pseudoresistor"
        tag5="Cadence w pseudoresistor"
        ini.defaultColors      = [color1,color3,color5]
    if (CHAPTER==42 or CHAPTER==44 or CHAPTER==51 or CHAPTER==52):
        for key in keys:
            if Trace1sel[key].label in selectionList:
                TracePlot[key]=Trace1sel[key]
                TracePlot[key].label=TracePlot[key].label.replace("Prelayout. Ideal resistor.",tag2)
                TracePlot[key].label=TracePlot[key].label.replace("Postlayout. Ideal resistor.",tag3)
                TracePlot[key].label=TracePlot[key].label.replace("Prelayout. Pseudoresistor.",tag4)
                TracePlot[key].label=TracePlot[key].label.replace("Postlayout. Whole circuit.",tag5)
    else:
        TracePlot=Trace1sel
        ini.defaultColors      = [color1,color2,color4,color3,'gray10','gray10',color5]
    SLiCAPPlotSweepPhase.axes[0][0].traces[0].label=SLiCAPPlotSweepPhase.axes[0][0].traces[0].label.replace("loopgain",tag1)
    traces2fig(TracePlot,SLiCAPPlotSweepPhase)
    SLiCAPPlotSweepPhase.plot()
    fig2html(SLiCAPPlotSweepPhase, 600)
    del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel, TracePlot
    ini.defaultColors      = ['r','b','g','c','m','y','k']

if P_Res_Ideal_4vs8_ResponsesDCAC:
    ini.plotFontSize=18
    htmlPage("DC pseudoresistance 8 vs 4 blocks")
    Trace1 = csv2traces('pseudoresistanceIdeal8Blocks.csv')
    Trace2 = csv2traces('pseudoresistanceIdeal4Blocks.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['ExperimentCompounREq (Vc=0.897) Y','ExperimentCompounREq (Vc=0.273) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("ExperimentCompounREq (", " ").replace(") Y","")+' 8 blocks'
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['ExperimentCompounREq (Vc=0.246) Y','ExperimentCompounREq (Vc=0.75) Y']:
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("ExperimentCompounREq (", " ").replace(") Y","")+' 4 blocks'
    Trace1sel.update(Trace2sel) #merge dictionaries
    print(list(Trace1sel.keys()))
    Graph = plot('LaTEX_PseudoResistor_DC_4vs8', 'DC response: 4 vs 8 blocks, Vclamp=1.65', 'semilogy', Trace1sel, xName = '$V_{DA}$', xUnits = 'V', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [-2,2] , yLim = [], show = False)
    fig2html(Graph, 600)
    del Trace1, Trace2, Trace1sel, Trace2sel
    Trace1 = csv2traces('pseudoresistanceIdeal8BlockAC.csv')
    Trace2 = csv2traces('pseudoresistanceIdeal4BlockAC.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['ReqAC (Vc=0.897) Y','ReqAC (Vc=0.273) Y']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("ReqAC (", " ").replace(") Y","")+' 8 blocks'
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ['ReqAC (Vc=0.246) Y','ReqAC (Vc=0.75) Y']:
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("ReqAC (", " ").replace(") Y","")+' 4 blocks'
    Trace1sel.update(Trace2sel) #merge dictionaries
    print(list(Trace1sel.keys()))
    Graph = plot('pseudoresistance48AC', 'AC response: 4 vs 8 blocks, Vclamp=1.65', 'log', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)


if PseudoResistor_DC_RealvsIdeal:
    ini.plotFontSize=18
    htmlPage("DC pseudoresistance Gain 245M")
    Trace1 = csv2traces('pseudoresistanceIdeal8Blocks.csv') #Ideal Vc
    Trace2 = Cadencecsv2traces('pseudoresistanceIdealBias.csv')   #Ideal Ibias
    Trace3 = Cadencecsv2traces('ReqFinalDC_full.csv') #With DAC
    Trace4 = Cadencecsv2traces('ReqStandardDC_extension.csv') #Replacing the DAC by an ideal source (standard measurements)

    keys=list(Trace1.keys())
    Trace1sel={}
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ['Req_BothBias (Vclamp=1.65,A=0,Ic=7.993965e-08) Y']:
            print(key)
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("Req_BothBias (Vclamp=1.65,A=0,Ic=7.993965e-08) Y","Ideal cell")
    Trace1sel.update(Trace2sel)

    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['Req (Vclamp=1.65,A=0,Ic=7.296614e-08) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("Req (Vclamp=1.65,A=0,Ic=7.296614e-08) Y","Real Implementation")
    Trace1sel.update(Trace4sel)
    Graph = plot('LaTEX_PseudoResistor_DC_Real_vs_Ideal', 'Comparison DC response pseudo-resistor, Vclamp=1.65', 'semilogy', Trace1sel, xName = '$V_{DA}$', xUnits = 'V', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)

if PseudoResistor_AC_RealvsIdeal:
    ini.plotFontSize=18
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
            Trace1sel[key].label=Trace1sel[key].label.replace("ReqAC (Vc=0.45) Y","Ideal cell")
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['ReqAC (Vclamp=1.65,A=0,Ic=7.296614e-08) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace("ReqAC (Vclamp=1.65,A=0,Ic=7.296614e-08) Y","Real implementation")
    Trace1sel.update(Trace4sel)
    Graph = plot('LaTEX_PseudoResistor_AC_Real_vs_Ideal', 'Comparison AC response pseudo-resistor, Vclamp=1.65', 'log', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)
    #del Trace1, Trace2, Trace3, Trace4, Trace1sel, Trace2sel, Trace3sel, Trace4sel


if ACLevels:
    ini.plotFontSize=18
    htmlPage("Pseudo_ACLevels")
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
    Graph = plot('LaTEX_PseudoResistor_AC_Levels_Prelayout', 'Pseudo-resistance AC response, Vclamp=1.65 - Prelayout', 'log', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [1,1e9] , yLim = [], show = False)
    fig2html(Graph, 600)
    del Trace4, Trace4sel

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
    ini.defaultColors      = ['r','c','b','m','g','y','k']
    Graph = plot('LaTEX_PseudoResistor_AC_Levels_Postlayout', 'Pseudo-resistance AC response, Vclamp=1.65 - Postlayout', 'log', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'Pseudo-resistance', yUnits = '$\Omega$', xLim = [1,1e9] , yLim = [], show = False)
    fig2html(Graph, 600)
    ini.defaultColors      = ['r','b','g','c','m','y','k']
    del Trace4, Trace4sel

if PseudoresistorBiasing:
    ini.plotFontSize=18
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
    Graph = plot('LaTEX_PseudoResistor_WbLb', 'Sizing of Wb & Lb.', 'semilogx', Trace1sel, xName = 'Ic', xUnits = 'A', yName = 'Vc', yUnits = 'V', xLim = [1e-9,329e-9] , yLim = [], show = False)
    fig2html(Graph, 600)

if ResistanceSelection:
    ini.plotFontSize=18
    htmlPage("Rps triode vs resistor")
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
    Graph = plot('LaTEX_PseudoResistor_Rps', 'Loop Gain Phase: resistor vs triode parallel transistors', 'semilogx', Trace1sel, xName = 'Frequency', xUnits = 'Hz', yName = 'Loop Gain Phase', yUnits = '$deg$', xLim = [1e-2,10e10] , yLim = [-200,200], show = False)
    fig2html(Graph, 600)

if Transient:
    ini.plotFontSize=18
    htmlPage("Transient Gain 200MV A")
    ini.defaultColors      = ['#00A6D6']
    Trace1 = Cadencecsv2traces('MidgainTransientPrelayout3.csv',limiter=0.2)
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['/out (Vclamp_delta=1.5,tr=0.1) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Vclamp_delta=1.5,tr=0.1) Y","Prelayout w Pseudoresistor (Gain=200M, Vclamp_delta=1.5,tr=0.1)")
    Graph = plot('LaTEX_PrelayoutResistor_Transient200M', 'Output signal amplitude (Gain=200MV/A)', 'lin', Trace1sel, xName = 't', xUnits = 's', yName = 'Amplitude', yUnits = '$V$', xLim = [0, 199] , yLim = [], show = True, xScale='m')
    fig2html(Graph, 600)
    del Trace1, Trace1sel, Graph
    ini.defaultColors      = ['#c3312f','#00a390','g','c','m','y','k']
    Trace1 = Cadencecsv2traces('Summary_Im_reference.csv',limiter=0.2)
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['/I8/MINUS (Im_pk=5e-09,Vclamp_delta=1.5,tr=0.1) Y"']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace('/I8/MINUS (Im_pk=5e-09,Vclamp_delta=1.5,tr=0.1) Y"',"Membrane current. (Im_pk=5e-09)")
    Graph = plot('LaTEX_PrelayoutResistor_TransientIm', 'Membrane current stimuli', 'lin', Trace1sel, xName = 't', xUnits = 's', yName = 'Amplitude', yUnits = '$A$', xLim = [0, 199] , yLim = [], show = True, xScale='m', yScale='n')
    fig2html(Graph, 600)
    del Trace1, Trace1sel, Graph
    ini.defaultColors      = ['#99d28c','#c3312f','c','m','y','k']
    Trace1 = Cadencecsv2traces('Summary_Vclamp_reference.csv',limiter=0.2)
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['/inB (Im_pk=5e-09,Vclamp_delta=1.5,tr=0.1) Y"']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace('/inB (Im_pk=5e-09,Vclamp_delta=1.5,tr=0.1) Y"',"Voltage clamp. (Vclamp_delta=1.5,tr=0.1)")
    Graph = plot('LaTEX_PrelayoutResistor_TransientVclamp', 'Voltage Clamp', 'lin', Trace1sel, xName = 't', xUnits = 's', yName = 'Amplitude', yUnits = '$V$', xLim = [0, 199] , yLim = [], show = True, xScale='m')
    fig2html(Graph, 600)
    del Trace1, Trace1sel, Graph
    ini.defaultColors      = ['r','b','g','c','m','y','k']

if NoiseChapter51pre:
    ini.plotFontSize=18
    htmlPage("Noise spectrum chapter prelayout 5 1")
    #Trace1 = Cadencecsv2traces('NoiseController.csv')
    Trace1 = csv2traces('onoise_Aitor.csv')
    Trace2 = Cadencecsv2traces('NoiseSpectrumPrelayout3.csv')
    #Trace3 = Cadencecsv2traces('NoiseSpectrumPostlayout.csv')
    keys=list(Trace1.keys())
    print(str(keys))
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
         if key in ['onoise (WRps=1e-06) Y']:
             Trace1sel[key]=Trace1[key]
             Trace1sel[key].label=Trace1sel[key].label.replace("onoise (WRps=1e-06) Y","Controller - Gain 200M")
             print(str(Trace1sel[key].label))
    Trace2sel={}
    keys=list(Trace2.keys())
    for key in keys:
        if key in ["onoise (Vclamp=1.65,A=0,Idac=9e-08) Y"]:
            Trace2sel[(key+'x')]=Trace2[key]
            Trace2sel[(key+'x')].label=Trace2sel[(key+'x')].label.replace("onoise (Vclamp=1.65,A=0,Idac=9e-08) Y","Total noise, Prelayout. - Gain 200M")
            print(str(Trace2sel[(key+'x')].label))
    Trace1sel.update(Trace2sel)
    print(list(Trace1sel.keys()))
    ini.defaultColors      = ['#00a390','#c3312f','g','c','m','y','k']
    Graph = plot('NoiseSpectrum51', 'Output referred noise spectrum', 'log', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Noise', yUnits = '$V^2/Hz$', xLim = [10**-3,20000] , yLim = [], show = False)
    #del Trace1, Trace2, Trace3, Trace1sel, Trace2sel, Trace3sel
    fig2html(Graph, 600)
    ini.defaultColors      = ['#eb7245','#f1be3e','#00a390','c','m','y','k']
    Trace1 = csv2traces('NoiseRms3.csv')
    Trace2 = csv2traces('RMSNoiseController.csv')
    Trace3 = csv2traces('RMSNoisePostLayout.csv')
    keys=list(Trace1.keys())
    print(keys)
    Trace1sel={}
    for key in keys:
        if key in ['leafValue( Vrms_onoise "A" 1 ) vs leafValue( ReqAC_Prelayout "A" 1 ) Y','leafValue( Vrms_onoise "A" 0 ) vs leafValue( ReqAC_Prelayout "A" 0 ) Y\n']:
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace('leafValue( Vrms_onoise "A" 0 ) vs leafValue( ReqAC_Prelayout "A" 0 ) Y\n',"Total noise HR mode - Prelayout")
            Trace1sel[key].label=Trace1sel[key].label.replace('leafValue( Vrms_onoise "A" 1 ) vs leafValue( ReqAC_Prelayout "A" 1 ) Y',"Total noise LR mode - Prelayout")
    keys=list(Trace2.keys())
    Trace2sel={}
    for key in keys:
        if key in ['Vrms_onoise (Vclamp=1.65) Y']:
            Trace2sel[key]=Trace2[key]
            Trace2sel[key].label=Trace2sel[key].label.replace("Vrms_onoise (Vclamp=1.65) Y","Controller - Prelayout")
            print(str(Trace2sel[key].label))
    Trace1sel.update(Trace2sel)
    GraphNoiseRMS = plot('NoiseRMS51', 'Output referred noise - RMS (DC-10kHz)', 'log', Trace1sel, xName = 'Gain', xUnits = 'V/A', yName = 'onoise', yUnits = 'Vrms',xLim = [] , yLim = [], show = False)
    fig2html(GraphNoiseRMS, 600)
    del Trace2, Trace2sel
    ini.defaultColors      = ['r','b','g','c','m','y','k']

if NoiseChapter51post:
    ini.plotFontSize=18
    htmlPage("Noise comparison")
    ini.defaultColors      = ['#00A6D6','#0066a2','g','c','m','y','k']
    Trace1 = Cadencecsv2traces('NoiseController.csv')
    Trace2 = Cadencecsv2traces('NoiseSpectrumPrelayout3.csv')
    Trace3 = Cadencecsv2traces('NoiseSpectrumPostlayout.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
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
    Trace1sel={}
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
    Graph = plot('Noisecomparison', 'Output referred noise spectrum', 'log', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Noise', yUnits = '$V^2/Hz$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)

    htmlPage("NOISE RMS COMPARISON")
    ini.defaultColors      = ['#0066a2','#00b7d3','#eb7245','#f1be3e','g','c','m','y','k']
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
    Trace1sel={}
    keys=list(Trace3.keys())
    Trace3sel={}
    for key in keys:
        if key in ['leafValue( Vrms_onoise "A" 1 ) vs leafValue( Req "A" 1 ) Y\n','leafValue( Vrms_onoise "A" 0 ) vs leafValue( Req "A" 0 ) Y']:
            Trace3sel[key]=Trace3[key]
            Trace3sel[key].label=Trace3sel[key].label.replace('leafValue( Vrms_onoise "A" 1 ) vs leafValue( Req "A" 1 ) Y',"Total noise LR mode - Postlayout")
            Trace3sel[key].label=Trace3sel[key].label.replace('leafValue( Vrms_onoise "A" 0 ) vs leafValue( Req "A" 0 ) Y',"Total noise HR mode - Postlayout")
            print(str(Trace3sel[key].label))
    Trace1sel.update(Trace3sel)
    GraphNoiseRMS = plot('NoiseRMScomparison', 'Output referred noise - RMS (DC-10kHz)', 'log', Trace1sel, xName = 'Gain', xUnits = 'V/A', yName = 'onoise', yUnits = 'Vrms',xLim = [] , yLim = [], show = False)
    fig2html(GraphNoiseRMS, 600)
    ini.defaultColors      = ['r','b','g','c','m','y','k']


if gainAC:
    ini.plotFontSize=18
    htmlPage("AC TRANSFER")
    Trace1 = Cadencecsv2traces('GainControllerMag.csv') #Ideal resistor
    Trace2 = Cadencecsv2traces('GainControllerMagPrelayout.csv')
    Trace3 = Cadencecsv2traces('GainControllerMagPostlayout.csv')
    keys=list(Trace1.keys())
    Trace1sel={}
    for key in keys:
        if key in ['Gain dB20 (Vclamp=1.65,Rfb=2000000) Y','Gain dB20 (Vclamp=1.65,Rfb=2e+08) Y','Gain dB20 (Vclamp=1.65,Rfb=7e+08) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("Gain dB20 (Vclamp=1.65,Rfb=2000000) Y","Gain Vclamp=1.65 Rfb=2M - Prelayout")
            Trace1sel[key].label=Trace1sel[key].label.replace("Gain dB20 (Vclamp=1.65,Rfb=2e+08) Y","Gain Vclamp=1.65 Rfb=200M - Prelayout")
            Trace1sel[key].label=Trace1sel[key].label.replace("Gain dB20 (Vclamp=1.65,Rfb=7e+08) Y","Gain Vclamp=1.65 Rfb=700M - Prelayout")
    Graph = plot('actransfer', 'Gain - AC transfer', 'semilogx', Trace1sel, xName = 'f', xUnits = 'Hz', yName = 'Gain', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)

if servoCompare:
    ini.plotFontSize=18
    htmlPage("Servo s")
    Trace4 = Cadencecsv2traces('ServoMagPrePostLayoutComparisonWholeCircuit.csv') #Replacing the DAC by an ideal source (standard measurements)
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['Servo dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=3.29e-07) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=1e-09) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=1.813836e-08) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y',"Gain=20G, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=1.813836e-08) Y","Gain=900M, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y","Gain=200M, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=3.29e-07) Y","Gain=60M, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=1e-09) Y',"Gain=60M, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y',"Gain=2M, LR mode")
    Graph = plot('servo_postlayout', 'Servo function - Post layout, Vclamp=1.65', 'semilogx', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'Servo Magnitude', yUnits = '$dB$', xLim = [] , yLim = [-20,5], show = False)
    fig2html(Graph, 600)
    del Trace4, Trace4sel
    Trace4 = Cadencecsv2traces('ServoMagPrelayout3.csv') #Replacing the DAC by an ideal source (standard measurements)
    keys=list(Trace4.keys())
    Trace4sel={}
    for key in keys:
        if key in ['Servo dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=3.29e-07) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y','Servo dB20 (Vclamp=1.65,A=1,Idac=1e-09) Y','Servo dB20 (Vclamp=1.65,A=0,Idac=1.813836e-08) Y']:
            print(key)
            Trace4sel[key]=Trace4[key]
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=0,Idac=1e-09) Y',"Gain=20G, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=1.813836e-08) Y","Gain=900M, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=9e-08) Y","Gain=200M, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace("Servo dB20 (Vclamp=1.65,A=0,Idac=3.29e-07) Y","Gain=60M, HR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=1e-09) Y',"Gain=60M, LR mode")
            Trace4sel[key].label=Trace4sel[key].label.replace('Servo dB20 (Vclamp=1.65,A=1,Idac=3.29e-07) Y',"Gain=2M, LR mode")
    Graph = plot('servo_prelayout', 'Servo function - Pre layout, Vclamp=1.65', 'semilogx', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'Servo Magnitude', yUnits = '$dB$', xLim = [] , yLim = [-20,5], show = False)
    fig2html(Graph, 600)
    del Trace4, Trace4sel

if PSRRiPlot:
    ini.plotFontSize=18
    htmlPage("PSRRi")
    Trace4 = Cadencecsv2traces('PSRRiPreLayout.csv')
    Trace3 = Cadencecsv2traces('PSRRiPostLayout.csv')
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
    Graph = plot('PSRRi_prelayout', 'PSRRi - Pre layout', 'semilogx', Trace4sel, xName = 'f', xUnits = 'Hz', yName = 'PSRRi', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
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
    Graph = plot('PSRRi_postlayout', 'PSRRi - Post layout', 'semilogx', Trace3sel, xName = 'f', xUnits = 'Hz', yName = 'PSRRi', yUnits = '$dB$', xLim = [] , yLim = [], show = False)
    fig2html(Graph, 600)

if PSRRvPlot:
    ini.plotFontSize=18
    htmlPage("PSRRv")
    Trace4 = Cadencecsv2traces('PSRRvPreLayout.csv')
    Trace3 = Cadencecsv2traces('PSRRvPostLayout.csv')
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

if Driving:
    ini.plotFontSize=18
    htmlPage("TDriving requirements")
    ini.defaultColors      = ['#00A6D6']
    Trace1 = Cadencecsv2traces('OutputStage_Transients.csv',limiter=300)
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['/out (Rfb=2e+08,L=3.5e-07,W=3.5e-07) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace("/out (Rfb=2e+08,L=3.5e-07,W=3.5e-07) Y","$V_{out}$")
    Graph = plot('LaTEX_PrelayoutResistor_Driving_V', 'Transient voltage', 'lin', Trace1sel, xName = 't', xUnits = 's', yName = 'Amplitude', yUnits = '$V$', xLim = [0, 299] , yLim = [], show = False, xScale = 'u')
    fig2html(Graph, 600)
    del Trace1, Trace1sel, Graph
    ini.defaultColors      = ['#c3312f','#00a390','#f1be3e','#99d28c','m','y','k']
    Trace1 = Cadencecsv2traces('OutputStage_Transients.csv',limiter=300)
    keys=list(Trace1.keys())
    print(str(keys))
    Trace1sel={}
    for key in keys:
        if key in ['/M1/D (Rfb=2e+08,L=3.5e-07,W=3.5e-07) Y','/C0/PLUS (Rfb=2e+08,L=3.5e-07,W=3.5e-07) Y','/I9/MINUS (Rfb=2e+08,L=3.5e-07,W=3.5e-07) Y']:
            print(key)
            Trace1sel[key]=Trace1[key]
            Trace1sel[key].label=Trace1sel[key].label.replace('/M1/D (Rfb=2e+08,L=3.5e-07,W=3.5e-07) Y',"$I_d$")
            Trace1sel[key].label=Trace1sel[key].label.replace('/C0/PLUS (Rfb=2e+08,L=3.5e-07,W=3.5e-07) Y',"$I_{load}$")
            Trace1sel[key].label=Trace1sel[key].label.replace('/I9/MINUS (Rfb=2e+08,L=3.5e-07,W=3.5e-07) Y',"$I_{bias}$")
    Graph = plot('LaTEX_PrelayoutResistor_Driving_I', 'Transient current', 'lin', Trace1sel, xName = 't', xUnits = 's', yName = 'Amplitude', yUnits = '$A$', xLim = [0, 299] , yLim = [], show = False, xScale = 'u', yScale='u')
    fig2html(Graph, 600)
    del Trace1, Trace1sel, Graph
    ini.defaultColors      = ['r','b','g','c','m','y','k']
