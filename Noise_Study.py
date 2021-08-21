#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from SLiCAP import *
from sympy.plotting import plot3d
from scipy.integrate import quad
from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt

def AssignProcessParameters(self):
    """
    Sets simplified EKV parameters

    """
    #######???????????????????????????????
    #COPY HERE PROCESS PARAMETERS AVAILABLE IN THE WIKI
    #ELSE, THE SCRIPT USES DEFAULT CMOS18 EKV PARAMETERS FROM BINKLEY
    #######???????????????????????????????
    pass

def LoadBasicCircuit(Name):
    """
    Loads SLiCAP netlist and assigns EKV paramaters

    """
    head2html('Circuit data')     #Creates HTML section for the netlist analysis results
    img2html(Name+'.svg',600)     #Plots the schematic into the HTML section
    makeNetlist(Name+'.asc',Name) #Creates a netlist .cir from a LTSpice schematic .asc
    netlist2html(Name+'.cir')     #Shows the netlist in the HTML section
    i1=instruction()              #Creates object to analize
    i1.setCircuit(Name+'.cir')    #Defines the netlist to be analyzed
    print("Circuit set")
    AssignProcessParameters(i1)   #Assigns the EKV process paramters to the elements of the netlist
    # Assigns model specefications
    i1.defPar('IG', IG)
    i1.defPar('Ce', Ce)
    i1.defPar('Ree', Ree)
    i1.defPar('Rs', Rs)
    return i1

def PlotAllSpectrum(self, Name, sel_gain=700e6, sel_W=30e-6, sel_I=2e-6, sel_L=2e-6, NEF=1, output=True, input=True, membrane=True):
    self.defPar('L', sel_L)      # Variable
    self.defPar('W', sel_W)      # Variable
    self.defPar('ID', sel_I)     # Variable
    self.defPar('gain', sel_gain)# Variable
    self.setSimType('numeric')   # Avoids unnecesary computation effort by replacing
                                 #with numeric values all variables whose value we know.
                                 #Applicable to all variables but W, ID and gain.
    self.setGainType('vi')       #to obtain transfer functions
    self.setDataType('noise')    #to calculate noise equations from noise sources at the .cir
    self.setSource('I1_XU1')     #the only noise source in the model is the controller's noise equivalent
    self.setDetector('V_out')    #the detector where we measure the output referred noise (input of the ADC = output of the amplifier)
    result = self.execute()      #obtain equations
    result.onoise=result.onoise*NEF #correct the output referred noise taking into account the NEF from the input stage biasing estimated in the script "Stages study"
    result.inoise=result.inoise*NEF #correct the input referred noise taking into account the NEF from the input stage biasing estimated in the script "Stages study"
    if output: #plot output referred spectrum
        figOnoiseSpec=plotSweep('onoise_spect'+Name+str(sel_gain),'Output referred noise sprectrum ('+Name+')', result, 8e-3, 2e4, 200, funcType='onoise', show=False, yLim = [10**-20,2*10**-9]) #plot the spectrum
        fig2html(figOnoiseSpec, 800)
    if input:  #plot input referred spectrum
        figInoiseSpec=plotSweep('inoise_spect'+Name+str(sel_gain),'Input referred noise sprectrum('+Name+')', result, 8e-3, 2e4, 200, funcType='inoise', show=False) #plot the spectrum
        fig2html(figInoiseSpec, 800)
    if membrane: #plot membrane referred spectrum
        self.setSource('I1')
        result2 = self.execute()
        result2.onoise=result.onoise*NEF
        figMnoiseSpec=plotSweep('mnoise_spect'+Name+str(sel_gain),'Membrane referred noise sprectrum('+Name+')', result2, 8e-3, 2e4, 200, funcType='inoise', show=False) #plot the spectrum
        fig2html(figMnoiseSpec, 800)
    return self

t1 = time()
prj = initProject('MSc Thesis')
print("Project inititated")
ini.MaximaTimeOut = 600         #Some extra time to allow the computer to calculate integrals.
ini.plotFontSize=14

#VARIABLES TO ADJUST DURING THIS SCRIPT:
Wopt=30e-6                      #Transistor width -> Obtained in Simulation 3
Iopt=2e-6                       #Drain current bias first stage -> Obtained in Simulation 4
L = 2e-6                        #Transistor length -> Obtained in Simulation 5

#VARIABLE OBTAINED IN SCRIPT "STAGES STUDY"
NEF=1.25

#SPECIFICATIONS
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


#EXECUTION PARAMETERS:
Name='NoiseController'          #Netlist selection (NoiseControler vs. Noise)
Optimization3D=False            #Generate 3-D graphs to pick Wopt
IDpick=True                     #Generate 2-D graphs to pick Iopt
lengthChecks=False              #Generate comparisons to pick L
NoiseReport=False               #Generate noise results summary (spectrums, rms)
if NoiseReport:
    Optimization3D=False
    IDpick=False
    lengthChecks=False
NoiseAnalysis=False             #Generate noise spectrum decomposition into colors

#SIMULATION #3
###############################################################################
###############################################################################
# Simulation considering only the controller noise source for its optimization.
# The model is simplified: the source impedance is formed by Ze and Rs.
# The gain is considered as a numeric variable
# The simulation is repeated for the two extreme gains and an intermediate one
# to check that the conclusions are not valid for just one gain.
###############################################################################
###############################################################################

#Establishes noise target based on ADC quantization noise
LSB=vRange/(2**res)                                             #LSB of the ADC
target=LSB/(sp.sqrt(12))*(sp.sqrt(noise_budget_vs_feedback))    #quantization noise * budget ratio by the controller
print("Qnoise: "+str(sp.N(LSB/(sp.sqrt(12)),4)))
print("Target controller: "+str(sp.N(target),4))

#Name='NoiseControllerPMOS'  #Activate to simulate for PMOS
htmlPage('Controller noise optimization')
i1=LoadBasicCircuit(Name)                                       #Loads the noise model
i1.defPar('L', L)
ID  = sp.Symbol('ID')                                           #Creates the symbolic variables
W  = sp.Symbol('W')
elementData2html(i1.circuit)
params2html(i1.circuit)
print("Parameters assigned")
##################

#Setting simulation config
print("Calculating noise")
i1.setSimType('numeric')   #Avoids unnecesary computation effort by replacing
                           #with numeric values all variables whose value we know.
                           #Applicable to all variables but W, ID and gain.
i1.setGainType('vi')       #Obtains transfer functions
i1.setDataType('noise')    #Calculates noise equations from noise sources at the .cir
#i1.setSource('I1_XU1')    #Not necessary: output-referred calculation
i1.setDetector('V_out')    #The detector where we measure the output referred noise (input of the ADC = output of the amplifier)


if Optimization3D:

    #################### SIMULATION 3a: Gain 20GV/A################
    ################################################################

    print('Evaluating for gain=20GV/A')
    head2html("Graph for gain=20GV/A")

    i1.defPar('gain', max_gain)
    t1 = time()
    result = i1.execute()   #Calculate all noise expressions for gain=20GV/A
    t2 = time()
    print(t2-t1, "s\n")
    t3 = time()


    """
    #SYMBOLIC APPROACH (WORKING IN PREVIOUS SLICAP VERSIONS)
    print("Simplifying expressions")
    # Simplify all the noise expressions by assuming positive parameters
    # After simplification remove the assumptions
    result.onoise = sp.simplify(assumePosParams(result.onoise))
    result.onoise = sp.N(clearAssumptions(result.onoise))
    result.onoise = sp.N(sp.simplify(result.onoise))
    print("simplifying results")
    for key in list(result.onoiseTerms.keys()):
        result.onoiseTerms[key] = sp.simplify(assumePosParams(result.onoiseTerms[key]))
        result.onoiseTerms[key] = sp.N(clearAssumptions(result.onoiseTerms[key]))
        print(key, result.onoiseTerms[key], '\n')

    RMSnoise = rmsNoise(result, 'onoise', f_min, f_max) #Calculate RMS value from spectrum noise result
    print("Noise converted to Vrms at the output")

    head2html('RMS Noise in terms of W and ID')
    expr2html(RMSnoise,"V^2/Hz") #Show results in HTML
    print(RMSnoise)

    head2html("Graph for gain=20G") #Plot result in 3D
    A=plot3d(RMSnoise, (W,180e-9,10000e-6), (ID,1e-6,0.01), title="$Vonoise_{rms}$ [$V^2$/Hz]", xlabel="W [m]", ylabel="ID [A]", cmap="RdYlBu_r", show=False)
    B=plot3d(RMSnoise, (W,20e-6,60e-6), (ID,50e-6,600e-6), title="$Vonoise_{rms}$ [$V^2$/Hz]", xlabel="W [m]", ylabel="ID [A]", cmap="RdYlBu_r", show=False) #Zoom version

    #Show plots in HTML
    A.save('./img/fullnoise_simplifiedmodel_maxgain.svg')
    img2html('fullnoise_simplifiedmodel_maxgain.svg',600)
    B.save('./img/fullnoise_simplifiedmodel_maxgain_zoom.svg')
    img2html('fullnoise_simplifiedmodel_maxgain_zoom.svg',600)
    """

    print("Creating plot")  #Plot the 3-D graph of onoise rms vs ID vs W, gain=200GV/A
    numW = 100
    numI = 10
    widths = np.linspace(1e-6, 100e-6, num=numW) # X variable
    IDs    = np.linspace(1e-6, 10e-6, num=numI) # Y variable
    X, Y   = np.meshgrid(widths, IDs)
    X      = X*1e6 # um
    Y      = Y*1e6 # uA
    Z      = np.zeros([numI, numW])
    # Create a function for numeric substitution of grid point values
    onoise = sp.lambdify([ID, W], result.onoise)
    for i in range(numI):
        for j in range(numW):
            # Convert sympy function to numpy function for current grid point
            onoiseF = sp.lambdify(ini.frequency, onoise(IDs[i], widths[j]))
            # Calculate RMS value for current grid point
            Z[i, j] = np.sqrt(quad(onoiseF, f_min, f_max)[0])*1e3*NEF # mV
    fig = plt.figure(figsize=(15,10))
    ax =  fig.add_subplot(1,1,1, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis', label='RMS noise [mV]')
    ax.set_xlabel('$W$ [um]', fontsize=18, rotation=150)
    ax.set_ylabel('$ID$ [uA]', fontsize=18)
    ax.set_zlabel('$V_{onoise}$ [mV(rms)]', fontsize=18, rotation=60)
    plt.show()
    t4 = time()
    print(t4-t3, "s\n")

    fig.savefig('./img/fullnoise_simplifiedmodel_maxgain.svg')
    img2html('fullnoise_simplifiedmodel_maxgain.svg',600)

    #################### SIMULATION 3b: Gain 200MV/A################
    ################################################################

    print('Evaluating for gain=200MV/A')
    head2html("Graph for gain=200MV/A")
    i1.defPar('gain', m_gain)
    t1 = time()
    result = i1.execute()   #Calculate all noise expressions for gain=200MV/A
    t2 = time()
    print(t2-t1, "s\n")
    t3 = time()

    print("Creating plot")  #Plot the 3-D graph of onoise rms vs ID vs W, gain=200MV/A
    numW = 100
    numI = 10
    widths = np.linspace(1e-6, 100e-6, num=numW) # X variable
    IDs    = np.linspace(1e-6, 10e-6, num=numI) # Y variable
    X, Y   = np.meshgrid(widths, IDs)
    X      = X*1e6 # um
    Y      = Y*1e6 # mA
    Z      = np.zeros([numI, numW])
    # Create a function for numeric substitution of grid point values
    onoise = sp.lambdify([ID, W], result.onoise)
    for i in range(numI):
        for j in range(numW):
            # Convert sympy function to numpy function for current grid point
            onoiseF = sp.lambdify(ini.frequency, onoise(IDs[i], widths[j]))
            # Calculate RMS value for current grid point
            Z[i, j] = np.sqrt(quad(onoiseF, f_min, f_max)[0])*1e3*NEF # mV
    fig = plt.figure(figsize=(15,10))
    ax =  fig.add_subplot(1,1,1, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis', label='RMS noise [mV]')
    ax.set_xlabel('$W$ [um]', fontsize=18, rotation=150)
    ax.set_ylabel('$ID$ [uA]', fontsize=18)
    ax.set_zlabel('$V_{onoise}$ [mV(rms)]', fontsize=18, rotation=60)
    plt.show()
    t4 = time()
    print(t4-t3, "s\n")

    fig.savefig('./img/fullnoise_simplifiedmodel_mgain.svg')
    img2html('fullnoise_simplifiedmodel_mgain.svg',600)

    #################### SIMULATION 3c: Gain 2MV/A##################
    ################################################################

    print('Evaluating for gain=2MV/A')
    head2html("Graph for gain=2MV/A")
    i1.defPar('gain', min_gain)
    t1 = time()
    result = i1.execute()   #Calculate all noise expressions for gain=2MV/A
    t2 = time()
    print(t2-t1, "s\n")
    t3 = time()

    print("Creating plot") #Plot the 3-D graph of onoise rms vs ID vs W, gain=2MV/A
    numW = 100
    numI = 10
    widths = np.linspace(1e-6, 100e-6, num=numW) # X variable
    IDs    = np.linspace(1e-6, 10e-6, num=numI) # Y variable
    X, Y   = np.meshgrid(widths, IDs)
    X      = X*1e6 # um
    Y      = Y*1e6 # uA
    Z      = np.zeros([numI, numW])
    # Create a function for numeric substitution of grid point values
    onoise = sp.lambdify([ID, W], result.onoise)
    for i in range(numI):
        for j in range(numW):
            # Convert sympy function to numpy function for current grid point
            onoiseF = sp.lambdify(ini.frequency, onoise(IDs[i], widths[j]))
            # Calculate RMS value for current grid point
            Z[i, j] = np.sqrt(quad(onoiseF, f_min, f_max)[0])*1e3*NEF # mV
    fig = plt.figure(figsize=(15,10))
    ax =  fig.add_subplot(1,1,1, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis', label='RMS noise [mV]')
    ax.set_xlabel('$W$ [um]', fontsize=18, rotation=150)
    ax.set_ylabel('$ID$ [uA]', fontsize=18)
    ax.set_zlabel('$V_{onoise}$ [mV(rms)]', fontsize=18, rotation=60)
    plt.show()
    t4 = time()
    print(t4-t3, "s\n")

    fig.savefig('./img/fullnoise_simplifiedmodel_mingain.svg')
    img2html('fullnoise_simplifiedmodel_mingain.svg',600)




################### ASSIGN Wopt ACCORDINGLY #############################

if IDpick:

    #SIMULATION #4
    ###############################################################################
    ###############################################################################
    # Simulation considering only the controller noise source.
    # The model is simplified: the source impedance is formed by Ze and Rs.
    # The gain is the worst case: maximum gain. (In the second round a quick check for intermediate and minimum gain are also included)
    # W is set to an optimum value based on the results of Simulation #3
    # ID is sweep until meeting the noise specs or considering a trade-off between energy expenditure and noise improvement.
    ###############################################################################
    ###############################################################################


    head2html("RMS Noise with Wopt, worst case gain (20GV/A)")

    i1.defPar('W', Wopt)       #Optimum with obtained in Simulation #3
    i1.defPar('gain', max_gain)
    ID  = sp.Symbol('ID')      #Parameter to sweep

    #Setting simulation config
    i1.setSimType('numeric')   #Avoids unnecesary computation effort by replacing
                               #with numeric values all variables whose value we know.
                               #Applicable to all variables but W, ID and gain.
    i1.setGainType('vi')       #Obtains transfer functions
    i1.setDataType('noise')    #Calculates noise equations from noise sources at the .cir
    #i1.setSource('I1_XU1')    #Not necessary: output-referred calculation
    i1.setDetector('V_out')    #The detector where we measure the output referred noise (input of the ADC = output of the amplifier)

    result = i1.execute()

    print("Creating plot")
    numI  = 100
    IDS   = np.linspace(0.1e-6, 10e-6, num=numI) # Y variable
    X     = IDS*1e3 # mA
    Y     = np.zeros([numI])
    # Create a function for numeric substitution of grid point values
    onoise = sp.lambdify([ID], result.onoise)
    for i in range(numI):
            # Convert sympy function to numpy function for current grid point
            onoiseF = sp.lambdify(ini.frequency, onoise(IDS[i]))
            # Calculate RMS value for current grid point
            Y[i] = np.sqrt(quad(onoiseF, f_min, f_max)[0])*1e3*NEF # mV
    fig = plt.figure(figsize=(10,10))
    ax =  fig.add_subplot(1,1,1)
    ax.plot(X, Y, label='RMS noise [mV]')
    ax.title.set_text('Gain='+str(i1.getParValue('gain')))
    ax.set_xlabel('$ID$ [mA]', fontsize=12)
    ax.set_ylabel('$V_{onoise}$ [mV(rms)]', fontsize=12)
    plt.grid(which='both', axis='both')
    plt.show()

    fig.savefig('./img/setIDplot.svg')
    img2html('setIDplot.svg',600)

    """
    #SYMBOLIC APPROACH (WORKING IN PREVIOUS SLICAP VERSIONS)
    # Simplify all the noise expressions by assuming positive parameters
    # After simplification remove the assumptions
    result.onoise = sp.simplify(assumePosParams(result.onoise))
    result.onoise = clearAssumptions(result.onoise)
    for key in list(result.onoiseTerms.keys()):
        result.onoiseTerms[key] = sp.simplify(assumePosParams(result.onoiseTerms[key]))
        result.onoiseTerms[key] = clearAssumptions(result.onoiseTerms[key])
    RMSnoise = rmsNoise(result, 'onoise', f_min, f_max) #Calculate RMS value from spectrum noise result
    print(RMSnoise)
    setIDplot=sp.plot(target,RMSnoise,(ID,5e-1,5e-8), xscale ='log', yscale ='linear',xlabel='ID[A]', ylabel='onoise_RMS [Vrms]', legend=False, show='False')
    setIDplot[0].line_color = 'red'
    setIDplot[0].label = 'target controller noise'
    setIDplot[1].line_color = 'blue'
    setIDplot[1].label = 'estimated controller noise (CS stage)'
    setIDplot.legend=True
    setIDplot.save('./img/setIDplot.svg')
    img2html('setIDplot.svg',600)
    """

    head2html("RMS Noise with Wopt, intermediate case gain (200MV/A)")

    i1.defPar('W', Wopt)       #Optimum with obtained in Simulation #3
    i1.defPar('gain', m_gain)
    ID  = sp.Symbol('ID')      #Parameter to sweep

    #Setting simulation config
    i1.setSimType('numeric')   #Avoids unnecesary computation effort by replacing
                               #with numeric values all variables whose value we know.
                               #Applicable to all variables but W, ID and gain.
    i1.setGainType('vi')       #Obtains transfer functions
    i1.setDataType('noise')    #Calculates noise equations from noise sources at the .cir
    #i1.setSource('I1_XU1')    #Not necessary: output-referred calculation
    i1.setDetector('V_out')    #The detector where we measure the output referred noise (input of the ADC = output of the amplifier)

    result = i1.execute()

    print("Creating plot")
    numI  = 100
    IDS     = np.linspace(0.1e-6, 10e-6, num=numI) # Y variable
    X     = IDS*1e3 # mA
    Y     = np.zeros([numI])
    # Create a function for numeric substitution of grid point values
    onoise = sp.lambdify([ID], result.onoise)
    for i in range(numI):
            # Convert sympy function to numpy function for current grid point
            onoiseF = sp.lambdify(ini.frequency, onoise(IDS[i]))
            # Calculate RMS value for current grid point
            Y[i] = np.sqrt(quad(onoiseF, f_min, f_max)[0])*1e3*NEF # mV
    fig = plt.figure(figsize=(10,10))
    ax =  fig.add_subplot(1,1,1)
    ax.plot(X, Y, label='RMS noise [mV]')
    ax.title.set_text('Gain='+str(i1.getParValue('gain')))
    ax.set_xlabel('$ID$ [mA]', fontsize=12)
    ax.set_ylabel('$V_{onoise}$ [mV(rms)]', fontsize=12)
    plt.grid(which='both', axis='both')
    plt.show()

    fig.savefig('./img/setIDplotmgain.svg')
    img2html('setIDplotmgain.svg',600)

    head2html("RMS Noise with Wopt, minimum (2MV/A)")

    i1.defPar('W', Wopt)        #Optimum with obtained in Simulation #3
    i1.defPar('gain', min_gain)
    ID  = sp.Symbol('ID')       #Parameter to sweep

    #Setting simulation config
    i1.setSimType('numeric')   #To avoid unnecesary computation effort by replacing
                               #with numeric values all variables whose value we know.
                               #Applicable to all variables but W, ID and gain.
    i1.setGainType('vi')       #to obtain transfer functions
    i1.setDataType('noise')    #to calculate noise equations from noise sources at the .cir
    #i1.setSource('I1_XU1')     #the only noise source in the model is the controller's noise equivalent
    i1.setDetector('V_out')    #the detector where we measure the noise (input of the ADC = output of the amplifier)

    result = i1.execute()

    print("Creating plot")
    numI  = 100
    IDS     = np.linspace(0.1e-6, 10e-6, num=numI) # Y variable
    X     = IDS*1e3 # mA
    Y     = np.zeros([numI])
    # Create a function for numeric substitution of grid point values
    onoise = sp.lambdify([ID], result.onoise)
    for i in range(numI):
            # Convert sympy function to numpy function for current grid point
            onoiseF = sp.lambdify(ini.frequency, onoise(IDS[i]))
            # Calculate RMS value for current grid point
            Y[i] = np.sqrt(quad(onoiseF, f_min, f_max)[0])*1e3*NEF # mV
    fig = plt.figure(figsize=(10,10))
    ax =  fig.add_subplot(1,1,1)
    ax.plot(X, Y, label='RMS noise [mV]')
    ax.title.set_text('Gain='+str(i1.getParValue('gain')))
    ax.set_xlabel('$ID$ [mA]', fontsize=12)
    ax.set_ylabel('$V_{onoise}$ [mV(rms)]', fontsize=12)
    plt.grid(which='both', axis='both')
    plt.show()

    fig.savefig('./img/setIDplotmingain.svg')
    img2html('setIDplotmingain.svg',600)

    ################### ASSIGN Iopt ACCORDINGLY #############################

if lengthChecks:

    #SIMULATION #5 Plotting the noise spectrum for different lengths (and adapting its Wopt accordingly)
    ###############################################################################
    ###############################################################################
    # Simulation considering the controller as the only noise source.
    # The model is simplified: the source impedance is formed by Ze and Rs.
    # The gain is numeric (worst case, max_gain).
    # W is set to an optimum value based on the results of Simulation #3 for different lengths.
    # ID is set to 10^-5 since it is a good trade-off between noise reduction and efficiency.
    ###############################################################################
    ###############################################################################
    i1.defPar('ID', Iopt)           #Optimum with obtained in Simulation #4
    i1.defPar('gain', max_gain)
    print("Parameters assigned")
    ini.stepFunction = False        #Fill the matrix with stepped numeric values

    #THIS TWO VECTORS SAVE THE Wopt AND IDopt FOR EXECUTIONS WITH DIFFERENT L TO COMPARE ALL OF THEM AND DECIDE UPON Lopt
    op_lengths=[350e-9,1e-6,2e-6,3.5e-6,5e-6]
    op_widths=[200e-6, 80e-6, 40e-6, 20e-6, 17.5e-6]

    stepArray=[op_widths, op_lengths]
    i1.setGainType('vi')       #generate laplace transfers
    i1.setDataType('noise')    #to calculate noise equations from noise sources at the .cir
    i1.setDetector('V_out')    #the detector where we measure the noise (input of the ADC = output of the amplifier)
    i1.setSimType('numeric')   #use numeric simulations
    i1.setStepMethod('array')  #try all the combinations
    i1.setStepVars(['W','L'])
    i1.setStepArray(stepArray)
    i1.stepOn()
    result=i1.execute()
    figOnoiseC=plotSweep('OnoiseC','ADC-referred noise sprectrum', result, f_min, f_max, 200, funcType='onoise', show='False') #plot the spectrums
    head2html('Noise spectrum (only controller)')
    fig2html(figOnoiseC, 800)
    stepArray2html(i1.stepVars, i1.stepArray)
    text2html('gain='+str(format(max_gain,'.1E'))+'V/A')
    text2html('ID='+str(format(Iopt,'.1E'))+'A')


################### NOISE REPORT AND COMPARISON WITH HAM ET AL ESTIMATED VALUE #############################

if NoiseReport:
    htmlPage("Report noise estimations")

    i1.defPar('W', Wopt)     #Optimum width obtained in Simulation #3
    i1.defPar('ID', Iopt)    #Optimum ID obtained in Simulation #4
    gain_sel  = sp.Symbol('gain')  #Parameter to sweep

    #Setting simulation config
    i1.setSimType('numeric')   #Avoids unnecesary computation effort by replacing
                               #with numeric values all variables whose value we know.
                               #Applicable to all variables but W, ID and gain.
    i1.setGainType('vi')       #Obtains transfer functions
    i1.setDataType('noise')    #Calculates noise equations from noise sources at the .cir
    i1.setSource('I1_XU1')     #the only noise source in the model is the controller's noise equivalent
    i1.setDetector('V_out')    #The detector where we measure the output referred noise (input of the ADC = output of the amplifier)
    resultA = i1.execute()

    i1.setSource('I1')
    resultB = i1.execute()

    ngtbc=100 #Number of gains to be checked
    gain_sel_lists=np.geomspace(min_gain, max_gain, num=ngtbc, endpoint=True).tolist()
    X=np.geomspace(min_gain/pow(10,6), max_gain/pow(10,6), num=ngtbc, endpoint=True).tolist()
    MyYo  = np.zeros([ngtbc])
    MyYi  = np.zeros([ngtbc])
    MyYm  = np.zeros([ngtbc])
    HamY  = np.zeros([ngtbc])
    HamYm  = np.zeros([ngtbc])

    # Create a function for numeric substitution of grid point values
    onoise = sp.lambdify([gain_sel], resultA.onoise)
    inoise = sp.lambdify([gain_sel], resultA.inoise)
    mnoise = sp.lambdify([gain_sel], resultB.inoise)


    for i in range(ngtbc):
            # Convert sympy function to numpy function for current grid point
            onoiseF = sp.lambdify(ini.frequency, onoise(gain_sel_lists[i]))
            inoiseF = sp.lambdify(ini.frequency, inoise(gain_sel_lists[i]))
            mnoiseF = sp.lambdify(ini.frequency, mnoise(gain_sel_lists[i]))
            # Calculate RMS value for current grid point
            MyYo[i] = np.sqrt(quad(onoiseF, f_min, f_max)[0])*1e3*NEF # mV
            #text2html('Output referred noise (DC-10kHz) [mVrms]: '+str(MyYo[i]))
            MyYi[i] = np.sqrt(quad(inoiseF, f_min, f_max)[0])*1e9*NEF # nV
            #text2html('Input referred noise (DC-10kHz) [nArms]: '+str(MyYi[i]))
            MyYm[i] = np.sqrt(quad(mnoiseF, f_min, f_max)[0])*1e12*NEF # pV
            #text2html('Membrane referred noise (DC-10kHz) [pArms]: '+str(MyYm[i]))
            HamY[i] = np.sqrt(quad(inoiseF, 1, 4700)[0])*1e12*NEF # pV
            #text2html('Input referred noise (1Hz-4.7kHz as measured by Ham): [pArms]'+str(HamY[i]))
            HamYm[i] = np.sqrt(quad(mnoiseF, 1, 4700)[0])*1e12*NEF # pV
            #text2html('Membrane referred noise (1Hz-4.7kHz as measured by Ham) [pArms]: '+str(HamYm[i]))

    head2html('Vrms noise vs gain (NEF included)')

    text2html('Output referred noise (DC-10kHz) [mVrms]') #PLOT RMS NOISE (OUTPUT REFERRED)
    fig = plt.figure(figsize=(15,10))
    ax =  fig.add_subplot(1,1,1)
    ax.plot(X, MyYo, label='Output referred noise (DC-10kHz) [mVrms]')
    ax.title.set_text('Output referred noise (DC-10kHz) [mVrms]')
    ax.set_xlabel('Gain [MV/A]', fontsize=12)
    ax.set_ylabel('$V_{onoise}$ [mV(rms)]', fontsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.grid(which='both', axis='both')
    plt.show()
    fig.savefig('./img/ONOISEREPORT.svg')
    img2html('ONOISEREPORT.svg',600)

    text2html('Input referred noise (DC-10kHz) [nArms]') #PLOT RMS NOISE (INPUT REFERRED)
    fig = plt.figure(figsize=(15,10))
    ax =  fig.add_subplot(1,1,1)
    ax.plot(X, MyYi, label='Input referred noise (DC-10kHz) [nArms]')
    ax.title.set_text('Input referred noise (DC-10kHz) [nArms]')
    ax.set_xlabel('Gain [MV/A]', fontsize=12)
    ax.set_ylabel('$V_{inoise}$ [nA(rms)]', fontsize=12)
    ax.set_xscale('log')
    plt.grid(which='both', axis='both')
    plt.show()
    fig.savefig('./img/INOISEREPORT.svg')
    img2html('INOISEREPORT.svg',600)

    text2html('Membrane referred noise (DC-10kHz) [pArms]') #PLOT RMS NOISE (MEMBRANE REFERRED)
    fig = plt.figure(figsize=(15,10))
    ax =  fig.add_subplot(1,1,1)
    ax.plot(X, MyYm, label='Membrane referred noise (DC-10kHz) [pArms]')
    ax.title.set_text('Membrane referred noise (DC-10kHz) [pArms]')
    ax.set_xlabel('Gain [MV/A]', fontsize=12)
    ax.set_ylabel('$V_{mnoise}$ [pA(rms)]', fontsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.grid(which='both', axis='both')
    plt.show()
    fig.savefig('./img/MNOISEREPORT.svg')
    img2html('MNOISEREPORT.svg',600)


    head2html('Spectrum plots for gain=700MV/A')  #PLOT OUTPUT REFERRED SPECTRUM (GAIN=700MV/A)
    i1.defPar('gain', 700e6)
    i1.setSource('I1_XU1')
    result = i1.execute()
    i1.setSource('I1')
    result2 = i1.execute()
    result.onoise=result.onoise*NEF
    result2.onoise=result.onoise*NEF
    figOnoiseSpec=plotSweep('onoise_spect','Output referred noise sprectrum. Gain=700M.', result, 8e-3, 2e4, 200, funcType='onoise', show=False) #plot the spectrum
    fig2html(figOnoiseSpec, 800)
    figInoiseSpec=plotSweep('inoise_spect','Input referred noise sprectrum. Gain=700M.', result, 8e-3, 2e4, 200, funcType='inoise', show=False) #plot the spectrum
    fig2html(figInoiseSpec, 800)
    figMnoiseSpec=plotSweep('mnoise_spect','Membrane referred noise sprectrum. Gain=700M.', result2, 8e-3, 2e4, 200, funcType='inoise', show=False) #plot the spectrum
    fig2html(figMnoiseSpec, 800)

    head2html('Spectrum plots for gain=200MV/A') #PLOT OUTPUT REFERRED SPECTRUMS (GAIN=200MV/A)
    i1.defPar('gain', 200e6)
    i1.setSource('I1_XU1')
    result = i1.execute()
    i1.setSource('I1')
    result2 = i1.execute()
    result.onoise=result.onoise*NEF
    result2.onoise=result.onoise*NEF
    ini.defaultColors=['k','#00b7d3','#99d28c','#f1be3e','#007188','k','#eb7245','#0066a2','#00a390','#82d7c6']
    figOnoiseSpec=plotSweep('onoise_spect_200','Output referred noise sprectrum. Gain=200M.', result, 8e-3, 2e4, 200, funcType='onoise', show=True) #plot the spectrum
    fig2html(figOnoiseSpec, 800)
    figInoiseSpec=plotSweep('inoise_spect_200','Input referred noise sprectrum. Gain=200M.', result, 8e-3, 2e4, 200, funcType='inoise', show=False) #plot the spectrum
    fig2html(figInoiseSpec, 800)
    figMnoiseSpec=plotSweep('mnoise_spect_200','Membrane referred noise sprectrum. Gain=200M.', result2, 8e-3, 2e4, 200, funcType='inoise', show=False) #plot the spectrum
    fig2html(figMnoiseSpec, 800)

    head2html('Spectrum plots for gain=2MV/A') #PLOT OUTPUT REFERRED SPECTRUMS (GAIN=2MV/A)
    i1.defPar('gain', 2e6)
    i1.setSource('I1_XU1')
    result = i1.execute()
    i1.setSource('I1')
    result2 = i1.execute()
    result.onoise=result.onoise*NEF
    ini.defaultColors=['k','#00b7d3','#99d28c','#f1be3e','#007188','k','#eb7245','#0066a2','#00a390','#82d7c6']
    figOnoiseSpec=plotSweep('onoise_spect_2','Output referred noise sprectrum. Gain=2M.', result, 8e-3, 2e4, 200, funcType='onoise', show=True) #plot the spectrum
    fig2html(figOnoiseSpec, 800)

    head2html('Spectrum plots for gain=20GV/A') #PLOT OUTPUT REFERRED SPECTRUMS (GAIN=20GV/A)
    i1.defPar('gain', 20e9)
    i1.setSource('I1_XU1')     #the only noise source in the model is the controller's noise equivalent
    result = i1.execute()
    i1.setSource('I1')
    result2 = i1.execute()
    result.onoise=result.onoise*NEF
    ini.defaultColors=['k','#00b7d3','#99d28c','#f1be3e','#007188','k','#eb7245','#0066a2','#00a390','#82d7c6']
    figOnoiseSpec=plotSweep('onoise_spect_20','Output referred noise sprectrum. Gain=20G.', result, 8e-3, 2e4, 200, funcType='onoise', show=True) #plot the spectrum
    fig2html(figOnoiseSpec, 800)


    head2html('Comparison Ham') #Compare results to Ham et all estimated values
    mnoiseF = sp.lambdify(ini.frequency,  result2.inoise)
    HamYm = np.sqrt(quad(mnoiseF, 1, 4700)[0])*1e12*NEF # pV
    text2html('Membrane referred noise (1Hz-4.7kHz as measured by Ham) [pArms]: '+str(HamYm))
    inoiseF = sp.lambdify(ini.frequency,  result.inoise)
    HamY = np.sqrt(quad(inoiseF, 1, 4700)[0])*1e12*NEF # pV
    text2html('Input referred noise (1Hz-4.7kHz as measured by Ham): [pArms]'+str(HamY))
    text2html("!! i don't know which one he measured")


if NoiseAnalysis: #Decompose the noise spectrum into
    Name="NoiseController_onlypink"  #Pink noise
    i1=LoadBasicCircuit(Name)
    Name="NoiseController_onlywhite" #White noise
    i2=LoadBasicCircuit(Name)
    Name="NoiseController_onlyblue"  #Blue noise
    i3=LoadBasicCircuit(Name)

    htmlPage("Gain 200MV A")
    Name="NoiseController_onlypink"
    ini.defaultColors=['#D33DC2']
    i1=PlotAllSpectrum(i1,Name,200e6, Wopt, Iopt, L, NEF, True, False, False)
    Name="NoiseController_onlywhite"
    ini.defaultColors=['#848484']
    i2=PlotAllSpectrum(i2,Name,200e6, Wopt, Iopt, L, NEF, True, False, False)
    Name="NoiseController_onlyblue"
    ini.defaultColors=['#00A6D6']
    i3=PlotAllSpectrum(i3,Name,200e6, Wopt, Iopt, L, NEF, True, False, False)

    htmlPage("Gain 2MV A")
    Name="NoiseController_onlypink"
    ini.defaultColors=['#D33DC2']
    i1=PlotAllSpectrum(i1,Name,2e6, Wopt, Iopt, L, NEF, True, False, False)
    Name="NoiseController_onlywhite"
    ini.defaultColors=['#848484']
    i2=PlotAllSpectrum(i2,Name,2e6, Wopt, Iopt, L, NEF, True, False, False)
    Name="NoiseController_onlyblue"
    ini.defaultColors=['#00A6D6']
    i3=PlotAllSpectrum(i3,Name,2e6, Wopt, Iopt, L, NEF, True, False, False)

    htmlPage("Gain 20GV A")
    Name="NoiseController_onlypink"
    ini.defaultColors=['#D33DC2']
    i1=PlotAllSpectrum(i1,Name,20e9, Wopt, Iopt, L, NEF, True, False, False)
    Name="NoiseController_onlywhite"
    ini.defaultColors=['#848484']
    i2=PlotAllSpectrum(i2,Name,20e9, Wopt, Iopt, L, NEF, True, False, False)
    Name="NoiseController_onlyblue"
    ini.defaultColors=['#00A6D6']
    i3=PlotAllSpectrum(i3,Name,20e9, Wopt, Iopt, L, NEF, True, False, False)
