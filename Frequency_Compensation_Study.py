#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from SLiCAP import *
from sympy.plotting import plot3d
from scipy.optimize import fminbound
from math import ceil, prod, log

def AssignProcessParameters(self):
    """
    Sets simplified EKV parameters

    """
    #######???????????????????????????????
    #COPY HERE PROCESS PARAMETERS AVAILABLE IN THE WIKI
    #ELSE, THE SCRIPT USES DEFAULT CMOS18 EKV PARAMETERS FROM BINKLEY
    #######???????????????????????????????
    pass

t1 = time()
prj = initProject('MSc Thesis')
print("Project inititated")
ini.MaximaTimeOut = 600         #Some extra time to allow the computer to calculate integrals.

# VARIABLES COMING FROM OTHER SCRIPTS
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

#Model parameters and specifications
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

#OTHER VARIABLES
ngtbc=3                         #Number of gains to be checked


#EXECUTION VARIABLES
circuits=["FreqComp_Both_NoBias"] #List of all circuits that will be checked. Options: FreqComp_Phantom_NoBias, FreqComp_PoleSplit_NoR_NoBias, FreqComp_PoleSplit_InclR_NoBias, FreqComp_Both
CpsDefaultBIG='5p'              #Pole splitting big capacitance  for gains of 200M or smaller
CpsDefaultSMALL='5p'            #Pole splitting small capacitance  for gains of 200M or larger
                                #Both Cps are equal because at the end a single capacitor is used
RpsDefault=200000               #Pole splitting resistance
Rph1Default=4000                #Phantom zero splitting resistance
PrintINFO_Before=False          #Print summary about stability before compensation
RootLocus_Feasibility=False     #Activate or deactivate large range root locus to determine feasibility
RootLocus_Value=True            #Activate or deactivate precise root locus to determine values of the compensation elements

for circuit in circuits:
    Name=circuit
    test=False
    makeNetlist(Name+'.asc',Name)       #create netlist
    i1=instruction()                    #create instruction object
    i1.setCircuit(Name+'.cir')          #load circuit
    print("-------------------")
    print("CIRCUIT SET: "+Name)
    AssignProcessParameters(i1)
    reference='Gm_M1_XU1'               #Sets the loop reference source for the asymtotic gain model
    i1.defPar('c_dg_XU1', 0)            #Eliminates the local feedback around the reference
    i1.defPar('c_dg_XU3', 0)            #Eliminates the local feedback around the reference
    i1.defPar('c_dg_XU2', 0)            #Eliminates the local feedback around the reference

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
        i1.defPar('range_ADC', rangeADC)
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

    # Assigns model specefications
    i1.defPar('IG', IG)
    i1.defPar('Ce', Ce)
    i1.defPar('Ree', Ree)
    i1.defPar('Rs', Rs)
    i1.defPar('Cload', Cload)
    print("Parameters assigned")

    gains=np.geomspace(min_gain, max_gain, num=ngtbc, endpoint=True).tolist() #generate list of gains to be checked

    #Print all circuit data
    htmlPage('Circuit data')            #Creates HTML section to show results
    head2html('Circuit diagram')
    img2html(Name + '.png', 1200)       #Plots the schematic into the HTML section
    print("Circuit data printed")
    netlist2html(Name + '.cir')         #Shows the netlist in the HTML section
    elementData2html(i1.circuit)        #print netlist elements and values
    params2html(i1.circuit)             #print execution parameters

    i1.setSource('I1')                  #the source is the membrane current source
    i1.setDetector('V_out')             #the detector is the output (input of the ADC)
    i1.setLGref(reference)              #reference source in the asymptotic model
    i1.setSimType('numeric')            #numeric simulation to substitute all variables with their values
    id=0

    #NO COMPENSATION IN THE BEGINNING
    try:
        i1.defPar('Cps', 0) #No pole splitting
        print('assigned value to Cps')
    except:
        print('no Cps')
    try:
        i1.defPar('Rps', 1e-6) #No pole splitting
        print('assigned value to Rps')
    except:
        print('no Rps')
    try:
        i1.defPar('Rph1', 1e-6) #No phantom zero
        print('assigned value to Rph1')
    except:
        print('no Rph1')

    for gain in gains:
        i1.defPar('gain', gain)
        id=id+1
        htmlPage('Compensation, gain: '+str(gain))
        if PrintINFO_Before:

            #Health check (the reuslts here needs to be equal to the ones in "Stages study" before compensation if cancelling the pole splitting and phantom zero were sucessful)
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

            #PRINT SUMMARY ABOUT STABILITY BEFORE COMEPNSATION (asymptotic gain model transfers Bode, bandwidth)
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

        if RootLocus_Feasibility:
            #ROOT LOCUS BRUTE-FORCE ANALYSIS: Freq. Compensation elements' values are swept in large ranges to determine if freq.comp. with that element is possible
            head2html('Compensation pole-splitting. Sweep.')
            gm4=i1.getParValue('g_m_XU4')
            try:
                if RpsDefault==None:
                    i1.defPar('Rps', 1/gm4)  #Assigns Rps (pole splitting resistance). Conventional value for the resistive element in the pole splitting: the inverse of the transconductance of the stage
                else:
                    i1.defPar('Rps', RpsDefault) #Assigns Rps (pole splitting resistance). Manually adjusted value to introduce a negative zero around the limit of the dominant frequencies
            except:
                print('no Rps')
            if gain<200e6:
                try:
                    i1.defPar('Cps', CpsDefaultBIG) #Assign (big) Cps (pole splitting capacitance) for gain<200e6
                except:
                    print('no Cps')
            else:
                try:
                    i1.defPar('Cps', CpsDefaultSMALL)  #Assign (small) Cps (pole splitting capacitance) for gain>=200e6
                except:
                    print('no Cps')
            try:
                i1.defPar('Rph1', Rph1Default) #Assign Rph (phantom zero resistor)
            except:
                print('no Rph1')
            i1.setGainType('gain')
            i1.setDataType('poles')
            gainp=i1.execute()
            i1.setStepMethod('log')
            #Draw root locus
            if (circuit=='FreqComp_Phantom_NoBias' or circuit=='FreqComp_Phantom'):
                i1.setStepVar('Rph1') #Sweeping Rph1
                i1.setStepStart(1e-6)
                i1.setStepStop(1e6)
                i1.setStepNum(100)
                i1.stepOn()
                figPolesServoPh1brute = plotPZ('FCpz'+str(id)+'brute'+circuit, 'poles of the gain', i1.execute(), show=False)
                figPolesServoPh1zbrute = plotPZ('FCpz'+str(id)+'brute'+'zoom'+circuit, 'poles of the gain (zoom)', i1.execute(), xmin=-1.2e5, xmax=1.2e5, ymin=-1.2e5, ymax=1.2e5, show=False)
                fig2html(figPolesServoPh1brute, 800)
                fig2html(figPolesServoPh1zbrute, 800)
            elif (circuit=='FreqComp_PoleSplit_NoR_NoBias' or circuit=='FreqComp_PoleSplit_NoR'):
                i1.setStepVar('Cps')  #Sweeping Cps
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
                i1.setStepVar('Rps') #Sweeping Rps
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
                i1.setStepVar('Rph1') #Sweeping Rph1
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

        if RootLocus_Value:
            #ROOT LOCUS FINE ADJUSTMENT: Freq. Compensation elements' values are swept with selected value to determine the exact value which is preferred
            head2html('Compensation pole-splitting. Selected values.')
            gm4=i1.getParValue('g_m_XU4')
            try:
                if RpsDefault==None:
                    i1.defPar('Rps', 1/gm4) #Assigns Rps (pole splitting resistance). Conventional value for the resistive element in the pole splitting: the inverse of the transconductance of the stage
                else:
                    i1.defPar('Rps', RpsDefault) #Assigns Rps (pole splitting resistance). Manually adjusted value to introduce a negative zero around the limit of the dominant frequencies
            except:
                print('no Rps')
            if gain<200e6:
                try:
                    i1.defPar('Cps', CpsDefaultBIG) #Assign (big) Cps (pole splitting capacitance) for gain<200e6
                except:
                    print('no Cps')
            else:
                try:
                    i1.defPar('Cps', CpsDefaultSMALL) #Assign (small) Cps (pole splitting capacitance) for gain>=200e6
                except:
                    print('no Cps')
            try:
                i1.defPar('Rph1', Rph1Default) #Assign Rph (phantom zero resistor)
            except:
                print('no Rph1')
            i1.setStepMethod('array')
            i1.setGainType('gain')
            if (circuit=='FreqComp_Phantom_NoBias' or circuit=='FreqComp_Phantom'):
                i1.setStepVars(['Rph1']) #Specific values of Rph1 to try
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
                i1.setStepVars(['Cps']) #Specific values of Cps to try
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
                i1.setStepVars(['Rps']) #Specific values of Rps to try
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
                i1.setStepVars(['Rph1']) #Specific values of Rph to try
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

        #PRINT SUMMARY ABOUT STABILITY AFTER COMEPNSATION (asymptotic gain model transfers Bode, bandwidth)
        gain=i1.getParValue('gain')
        try:
            if RpsDefault==None:
                i1.defPar('Rps', 1/gm4) #Assigns Rps (pole splitting resistance). Conventional value for the resistive element in the pole splitting: the inverse of the transconductance of the stage
            else:
                i1.defPar('Rps', RpsDefault) #Assigns Rps (pole splitting resistance). Manually adjusted value to introduce a negative zero around the limit of the dominant frequencies
        except:
            print('no Rps')
        if gain<200e6:
            try:
                i1.defPar('Cps', CpsDefaultBIG) #Assign (big) Cps (pole splitting capacitance) for gain<200e6
            except:
                print('no Cps')
        else:
            try:
                i1.defPar('Cps', CpsDefaultSMALL) #Assign (small) Cps (pole splitting capacitance) for gain>=200e6
            except:
                print('no Cps')
        try:
            i1.defPar('Rph1', Rph1Default)  #Assign Rph (phantom zero resistor)
        except:
            print('no Rph1')

        i1.setGainType('gain')
        i1.setDataType('poles')
        gainp=i1.execute()
        pz2html(gainp) #poles of the gain (stability check)
        #transfers:
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
        #Bode diagrams of the transfers:
        head2html('Bode plots')
        result = [asymptotic, gaing, loopgain, servo, direct]
        figdBmagPh = plotSweep('dBmagPhgain'+str(id)+circuit, 'dB magnitude', result, f_min, f_check, 100, funcType = 'dBmag', show=False)
        figPhasePh = plotSweep('phasePhgain'+str(id)+circuit,'Phase', result, f_min, f_check, 100, funcType = 'phase', show=False)
        fig2html(figdBmagPh, 800)
        fig2html(figPhasePh, 800)
        #Bandwidth
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
        #summary of the values used for compensation
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
