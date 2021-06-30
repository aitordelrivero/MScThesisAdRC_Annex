#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

SCRIPT #4: Design of frequency response.

"""


from SLiCAP import *
from sympy.plotting import plot3d
from scipy.optimize import fminbound
from math import ceil, prod, log

t1 = time()
prj = initProject('MSc Thesis')
print("Project inititated")
ini.MaximaTimeOut = 600         #Some extra time to allow the computer to calculate integrals.

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#List of all circuits that will be checked
circuits=["FreqComp_Both"] #"FreqComp_PoleSplit_InclR"
#One cps for gains of 200M or larger, and other for gains of 200M or lower. At the end I use only one value
CpsDefaultBIG='5p'
CpsDefaultSMALL='5p'
RpsDefault=200000 #None for 1/gm4
Rph1Default=4000
PrintPreCompensationInformation=False
BruteRootLocus=False
SelectiveRootLocus=True
#circuits=["FreqComp_PoleSplit_NoR_NoBias","FreqComp_PoleSplit_NoR","FreqComp_PoleSplit_InclR","FreqComp_PoleSplit_InclR_NoBias", "FreqComp_Phantom", "FreqComp_Phantom_NoBias","FreqComp_Both_NoBias"]
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# Transistor sizing (Instance i1)
Win=30e-6
IDin=2e-6
Lin = 2e-6
Wout=5e-6
Lout=350e-9
IDout=5e-6
Lcm=5e-6
Wcm=7.5e-6
Wb=5e-6
Lb=1.5e-6
Wbcm=8e-6
Lbcm=1.5e-6

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

for circuit in circuits:
    Name=circuit
    test=False

    makeNetlist(Name+'.asc',Name)       #create netlist
    i1=instruction()                    #create instruction object
    i1.setCircuit(Name+'.cir')          #load circuit
    print("-------------------")
    print("CIRCUIT SET: "+Name)

    #Loading model parameters into the .cir
    OhmNMOS=300
    OhmPMOS=300
    ###############????????????
    #Paste here process data at https://wiki-bsse.ethz.ch/display/DBSSEVLSI/EPFL-EKV+model+of+the+XFAB-xh018+process+for+simplified+design
    ###############????????????

    reference='Gm_M1_XU1'
    i1.defPar('c_dg_XU1', 0)
    i1.defPar('c_dg_XU3', 0)
    i1.defPar('c_dg_XU2', 0)

    try:
        i1.defPar('W', Win)
        i1.defPar('ID', IDin)
        i1.defPar('L', Lin)
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
        i1.defPar('IDas', IDin*2)
        i1.defPar('Was', Win*2)
        i1.defPar('Lin', Lin)
        print('input balanced CS parameters assigned')
    except:
        print('input antiseries not needed')
    try:
        i1.defPar('Wcm', Wcm)
        i1.defPar('Lcm', Lcm)
        i1.defPar('IDasPMOS', -IDin*2)
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
        i1.defPar('IDas', IDin*2)
        print('output bias equivalent parameters assigned')
    except:
        print('no input bias equivalent')

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
    elementData2html(i1.circuit)
    params2html(i1.circuit)

    i1.setSource('I1')       #the source is the membrane current source
    i1.setDetector('V_out')  #the detector is the output (input of the ADC)
    i1.setLGref(reference)
    i1.setSimType('numeric') #numeric simulation to substitute all variables with their values
    id=0

    try:
        i1.defPar('Cps', 0)
        print('assigned value to Cps')
    except:
        print('no Cps')
    try:
        i1.defPar('Rps', 1e-6)
        print('assigned value to Rps')
    except:
        print('no Rps')
    try:
        i1.defPar('Rph1', 1e-6)
        print('assigned value to Rph1')
    except:
        print('no Rph1')

    for gain in gains:
        i1.defPar('gain', gain)
        id=id+1
        htmlPage('Compensation, gain: '+str(gain))
        if PrintPreCompensationInformation:
            #PRINT SUMMARY ABOUT FREQUENCY RESPONSE BEFORE COMEPNSATION
            head2html('Gain poles before Ph')
            head2html('Check loopgain (must be equal to loopgain study depsite adding Chp1,Rph=0)')
            i1.setGainType('loopgain')
            i1.setDataType('pz')
            loopgainp=i1.execute()
            pz2html(loopgainp)
            i1.setGainType('gain')
            i1.setDataType('pz')
            gainpz=i1.execute()
            pz2html(gainpz)

            f_check=10e6
            i1.setGainType('gain')
            i1.setDataType('laplace')
            gaing = i1.execute()
            i1.setGainType('asymptotic')
            i1.setDataType('laplace')
            asymptotic = i1.execute()
            i1.setGainType('loopgain')
            i1.setDataType('laplace')
            loopgain = i1.execute()
            i1.setGainType('servo')
            i1.setDataType('laplace')
            servo = i1.execute()
            i1.setGainType('direct')
            i1.setDataType('laplace')
            direct = i1.execute()
            head2html('Bode plots')
            result = [asymptotic, gaing, loopgain, servo, direct]
            figdBmagPhDELETE = plotSweep('dBmagPhgainDELETE'+str(id)+circuit, 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
            figPhasePhDELETE = plotSweep('phasePhgainDELETE'+str(id)+circuit,'Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
            fig2html(figdBmagPhDELETE, 800)
            fig2html(figPhasePhDELETE, 800)
            print("Bode plots obtained")
            servoData = findServoBandwidth(loopgain.laplace)
            f_h = servoData['lpf']
            text2html('Bandwidth: '+str(format(f_h,'.3E'))+' Hz')
            del asymptotic, gaing, loopgain, servo, direct, gainpz, servoData

        if BruteRootLocus:
            #ROOT LOCUS BRUTE ANALYSIS
            head2html('Compensation pole-splitting. Sweep.')
            gm4=i1.getParValue('g_m_XU4')
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
            i1.setGainType('gain')
            i1.setDataType('poles')
            gainp=i1.execute()
            i1.setStepMethod('log')
            if (circuit=='FreqComp_Phantom_NoBias' or circuit=='FreqComp_Phantom'):
                i1.setStepVar('Rph1')
                i1.setStepStart(1e-6)
                i1.setStepStop(1e6)
                i1.setStepNum(100)
                i1.stepOn()
                figPolesServoPh1brute = plotPZ('FCpz'+str(id)+'brute'+circuit, 'poles of the gain', i1.execute(), show=False)
                figPolesServoPh1zbrute = plotPZ('FCpz'+str(id)+'brute'+'zoom'+circuit, 'poles of the gain (zoom)', i1.execute(), xmin=-1.2e5, xmax=1.2e5, ymin=-1.2e5, ymax=1.2e5, show=False)
                fig2html(figPolesServoPh1brute, 800)
                fig2html(figPolesServoPh1zbrute, 800)
            elif (circuit=='FreqComp_PoleSplit_NoR_NoBias' or circuit=='FreqComp_PoleSplit_NoR'):
                i1.setStepVar('Cps')
                i1.setStepStart(1e-15)
                i1.setStepStop(1e-9)
                i1.setStepNum(100)
                i1.stepOn()
                figPolesServoPh1brute = plotPZ('FCpz'+str(id)+'brute'+circuit, 'poles of the gain', i1.execute(), show=False)
                figPolesServoPh1zbrute = plotPZ('FCpz'+str(id)+'brute'+'zoom'+circuit, 'poles of the gain (zoom)', i1.execute(), xmin=-1.2e5, xmax=1.2e5, ymin=-1.2e5, ymax=1.2e5, show=False)
                fig2html(figPolesServoPh1brute, 800)
                fig2html(figPolesServoPh1zbrute, 800)
            elif (circuit=='FreqComp_PoleSplit_InclR_NoBias' or circuit=='FreqComp_PoleSplit_InclR'):
                i1.setDataType('zeros')
                i1.setStepVar('Rps')
                i1.setStepStart(1e-6)
                i1.setStepStop(2e6)
                i1.setStepNum(100)
                i1.stepOn()
                figPolesServoPh1brute = plotPZ('FCpz'+str(id)+'brute'+circuit, 'poles of the gain', i1.execute(), show=False)
                figPolesServoPh1zbrute = plotPZ('FCpz'+str(id)+'brute'+'zoom'+circuit, 'poles of the gain (zoom)', i1.execute(), xmin=-1.2e5, xmax=1.2e5, ymin=-1.2e5, ymax=1.2e5, show=False)
                fig2html(figPolesServoPh1brute, 800)
                fig2html(figPolesServoPh1zbrute, 800)
            elif (circuit=='FreqComp_Both_NoBias' or circuit=="FreqComp_Both"):
                i1.setDataType('zeros')
                i1.setStepVar('Rph1')
                i1.setStepStart(1e-6)
                i1.setStepStop(1e6)
                i1.setStepNum(100)
                i1.stepOn()
                figPolesServoPh1brute = plotPZ('FCpz'+str(id)+'brute'+circuit, 'poles of the gain', i1.execute(), show=False)
                figPolesServoPh1zbrute = plotPZ('FCpz'+str(id)+'brute'+'zoom'+circuit, 'poles of the gain (zoom)', i1.execute(), xmin=-1.2e5, xmax=1.2e5, ymin=-1.2e5, ymax=1.2e5, show=False)
                fig2html(figPolesServoPh1brute, 800)
                fig2html(figPolesServoPh1zbrute, 800)
            else:
                print('could not print root locus sweep')

        #SELECTED VALUES
        if SelectiveRootLocus:
            head2html('Compensation pole-splitting. Selected values.')
            gm4=i1.getParValue('g_m_XU4')
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
            i1.setStepMethod('array')
            i1.setGainType('gain')
            if (circuit=='FreqComp_Phantom_NoBias' or circuit=='FreqComp_Phantom'):
                i1.setStepVars(['Rph1'])
                stepArray=[1e-6, 1, 10, 100, 1000, 10000, 100000, 1e6, 1e7, 1e8]
                i1.setStepArray([stepArray])
                i1.stepOn()
                figPolesServoPh1 = plotPZ('FCpz'+str(id)+circuit, 'poles of the gain', i1.execute(), show=False)
                figPolesServoPh1z = plotPZ('polesServoPhZgain'+str(id)+'zoom'+circuit, 'poles of the gain (zoom)', i1.execute(), xmin=-3e5, xmax=3e5, ymin=-3e5, ymax=3e5, show=False)
                head2html('Gain poles after ph1')
                fig2html(figPolesServoPh1, 800)
                fig2html(figPolesServoPh1z, 800)
                stepArray2html(i1.stepVars, i1.stepArray)
                text2html('Rps='+str(format(i1.getParValue('Rps'),'.2E'))+'$\Omega$')
                i1.stepOff()
            elif (circuit=='FreqComp_PoleSplit_NoR_NoBias' or circuit=='FreqComp_PoleSplit_NoR'):
                i1.setStepVars(['Cps'])
                stepArray=[0, '1p', '2p', '3p', '4p', '5p', '6p']
                i1.setStepArray([stepArray])
                i1.stepOn()
                figPolesServoPh1 = plotPZ('FCpz'+str(id)+circuit, 'poles of the gain', i1.execute(), show=False)
                figPolesServoPh1z = plotPZ('polesServoPhZgain'+str(id)+'zoom'+circuit, 'poles of the gain (zoom)', i1.execute(), xmin=-1e3, xmax=1e3, ymin=-1e3, ymax=1e3, show=False)
                fig2html(figPolesServoPh1, 800)
                fig2html(figPolesServoPh1z, 800)
                stepArray2html(i1.stepVars, i1.stepArray)
                text2html('Rps='+str(format(i1.getParValue('Rps'),'.2E'))+'$\Omega$')
                i1.stepOff()
            elif (circuit=='FreqComp_PoleSplit_InclR_NoBias' or circuit=='FreqComp_PoleSplit_InclR'):
                i1.setDataType('zeros')
                i1.setStepVars(['Rps'])
                stepArray=[0.1, 1e3, 1e4, 1e5, 1e6, 1e7]
                i1.setStepArray([stepArray])
                i1.stepOn()
                figPolesServoPh1 = plotPZ('FCpz'+str(id)+circuit, 'poles of the gain', i1.execute(), xmin=-1e7, xmax=1e7, ymin=-1e7, ymax=1e7, show=False)
                figPolesServoPh1z = plotPZ('polesServoPhZgain'+str(id)+'zoom'+circuit, 'poles of the gain (zoom)', i1.execute(), xmin=-1e4, xmax=1e4, ymin=-1e4, ymax=1e4, show=False)
                fig2html(figPolesServoPh1, 800)
                fig2html(figPolesServoPh1z, 800)
                stepArray2html(i1.stepVars, i1.stepArray)
                i1.stepOff()
            elif (circuit=='FreqComp_Both_NoBias' or circuit=="FreqComp_Both"):
                i1.setDataType('zeros')
                i1.setStepVars(['Rph1'])
                stepArray=[0.1, 1e3, 1e4, 1e5, 1e6, 1e7]
                i1.setStepArray([stepArray])
                i1.stepOn()
                figPolesServoPh1 = plotPZ('FCpz'+str(id)+circuit, 'poles of the gain', i1.execute(), xmin=-1e7, xmax=1e7, ymin=-1e7, ymax=1e7, show=False)
                figPolesServoPh1z = plotPZ('polesServoPhZgain'+str(id)+'zoom'+circuit, 'poles of the gain (zoom)', i1.execute(), xmin=-1e4, xmax=1e4, ymin=-1e4, ymax=1e4, show=False)
                fig2html(figPolesServoPh1, 800)
                fig2html(figPolesServoPh1z, 800)
                stepArray2html(i1.stepVars, i1.stepArray)
                i1.stepOff()
            else:
                print('could not print root locus sweep')


        #FINAL SELECTION DIAGRAMS PLOT

        gain=i1.getParValue('gain')
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

        i1.setGainType('gain')
        i1.setDataType('poles')
        gainp=i1.execute()
        pz2html(gainp)

        f_check=10e6
        i1.setGainType('gain')
        i1.setDataType('laplace')
        gaing = i1.execute()
        i1.setGainType('asymptotic')
        i1.setDataType('laplace')
        asymptotic = i1.execute()
        i1.setGainType('loopgain')
        i1.setDataType('laplace')
        loopgain = i1.execute()
        i1.setGainType('servo')
        i1.setDataType('laplace')
        servo = i1.execute()
        i1.setGainType('direct')
        i1.setDataType('laplace')
        direct = i1.execute()
        head2html('Bode plots')
        result = [asymptotic, gaing, loopgain, servo, direct]
        figdBmagPh = plotSweep('dBmagPhgain'+str(id)+circuit, 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
        figPhasePh = plotSweep('phasePhgain'+str(id)+circuit,'Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
        fig2html(figdBmagPh, 800)
        fig2html(figPhasePh, 800)
        print("Bode plots obtained")
        servoData = findServoBandwidth(loopgain.laplace)
        f_h = servoData['lpf']
        text2html('Bandwidth: '+str(format(f_h,'.3E'))+' Hz')
        del asymptotic, gaing, loopgain, servo, direct, gainp, servoData

        #pole-zero study of the influence of compensation techniques on the loopgain
        head2html('Check loopgain after compensation')
        i1.setGainType('loopgain')
        i1.setDataType('pz')
        loopgainp=i1.execute()
        pz2html(loopgainp)
        try:
            text2html('Rps='+str(format(i1.getParValue('Rps'),'.2E'))+'$\Omega$, or none')
        except:
            print('no Rps')
        try:
            text2html('Cps='+str(format(i1.getParValue('Cps'),'.2E'))+'$\Omega$, or none')
        except:
            print('no Cps')
        try:
            text2html('Rph1='+str(format(i1.getParValue('Rph1'),'.2E'))+'$\Omega$, or none')
        except:
            print('no Rph1')
