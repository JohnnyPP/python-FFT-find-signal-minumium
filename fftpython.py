import numpy as np
import matplotlib.pyplot as plt
import math

phase=170
numberOfPoints=360
phaseRAD=np.deg2rad(phase)

#phase=90 gives for fft(sin(x+phaseRAD) fftPhaseHarmonic[1]=0
#the minimum is in 180 degrees

#phase=0 minimum is in 270 degrees fftPhaseHarmonic[1]=-90
#phase=10 minimum is in 260 degrees fftPhaseHarmonic[1]=-80
#phase=20 minimum is in 250 degrees fftPhaseHarmonic[1]=-70
#phase=30 minimum is in 240 degrees fftPhaseHarmonic[1]=-60
#phase=40 minimum is in 230 degrees fftPhaseHarmonic[1]=-50
#phase=50 minimum is in 220 degrees fftPhaseHarmonic[1]=-40
#phase=60 minimum is in 210 degrees fftPhaseHarmonic[1]=-30
#phase=70 minimum is in 200 degrees fftPhaseHarmonic[1]=-20
#phase=80 minimum is in 190 degrees fftPhaseHarmonic[1]=-10
#phase=90 minimum is in 180 degrees fftPhaseHarmonic[1]=0
#phase=100 minimum is in 170 degrees fftPhaseHarmonic[1]=10
#phase=110 minimum is in 160 degrees fftPhaseHarmonic[1]=20
#phase=120 minimum is in 150 degrees fftPhaseHarmonic[1]=30
#phase=130 minimum is in 140 degrees fftPhaseHarmonic[1]=40
#phase=140 minimum is in 130 degrees fftPhaseHarmonic[1]=50
#phase=150 minimum is in 120 degrees fftPhaseHarmonic[1]=60
#phase=160 minimum is in 110 degrees fftPhaseHarmonic[1]=70
#phase=170 minimum is in 100 degrees fftPhaseHarmonic[1]=80 
#phase=180 minimum is in 90 degrees fftPhaseHarmonic[1]=90
#phase=190 minimum is in 80 degrees fftPhaseHarmonic[1]=100
#phase=200 minimum is in 70 degrees fftPhaseHarmonic[1]=110
#phase=210 minimum is in 60 degrees fftPhaseHarmonic[1]=120
#phase=220 minimum is in 50 degrees fftPhaseHarmonic[1]=130
#phase=230 minimum is in 40 degrees fftPhaseHarmonic[1]=140
#phase=240 minimum is in 30 degrees fftPhaseHarmonic[1]=150
#phase=250 minimum is in 20 degrees fftPhaseHarmonic[1]=160
#phase=260 minimum is in 10 degrees fftPhaseHarmonic[1]=170
#phase=270 minimum is in 0 degrees fftPhaseHarmonic[1]=-180
#phase=280 minimum is in 350 degrees fftPhaseHarmonic[1]=-170
#phase=290 minimum is in 340 degrees fftPhaseHarmonic[1]=-160
#phase=300 minimum is in 330 degrees fftPhaseHarmonic[1]=-150
#phase=310 minimum is in 320 degrees fftPhaseHarmonic[1]=-140
#phase=320 minimum is in 310 degrees fftPhaseHarmonic[1]=-130
#phase=330 minimum is in 300 degrees fftPhaseHarmonic[1]=-120
#phase=340 minimum is in 290 degrees fftPhaseHarmonic[1]=-110
#phase=350 minimum is in 280 degrees fftPhaseHarmonic[1]=-100
#phase=360 minimum is in 270 degrees fftPhaseHarmonic[1]=-90

def FindMinimumPhase(FFT):
    """
    calculates the minimum phase for given angle obtained with
    the FFT algorithm. This calculations are taking into account sinus
    phase shift FFT(sin(x+phaseRAD)
    """    
    if FFT < 0 and FFT >= -180:
        if FFT != -180:
            minPhase=abs(FFT)+180
        else:
            minPhase=0
    
    if FFT >= 0 and FFT <= 180:
        minPhase=FFT-180
        minPhase=abs(minPhase)
        
    print minPhase
    return minPhase

x = np.linspace(0, 2*math.pi,numberOfPoints)   #create 360 angle points
xFFT = np.linspace(0,numberOfPoints/2+1,numberOfPoints/2+1)   
                                #the range is two times smaller +1 for RFFT
y = np.sin(x+phaseRAD)          #sinus signal without noise used for fit
ynoise = np.sin(x+phaseRAD)+0.5*np.random.normal(size=x.shape)  #noisy signal


dataLenghtFFT = len(ynoise)/2        #divide by 2 to satify rfft
                    # scale by the number of points so that
                    # the magnitude does not depend on the length 
                    # of the signal or on its sampling frequency  


calculatedFFT = np.fft.rfft(ynoise)     
#calculatedFFT = np.fft.rfft(y) 

#When the DFT is computed for purely real input, the output is 
#Hermite-symmetric, i.e. the negative frequency terms are just the complex 
#conjugates of the corresponding positive-frequency terms, and the 
#negative-frequency terms are therefore redundant. This function does not 
#compute the negative frequency terms, and the length of the transformed axis 
#of the output is therefore n//2+1

#When the input a is a time-domain signal and A = fft(a), np.abs(A) is its 
#amplitude spectrum and np.abs(A)**2 is its power spectrum. The phase spectrum 
#is obtained by np.angle(A).

amplitudeFFT = np.abs(calculatedFFT)    #calculates FFT amplitude from 
                                        #complex calculatedFFT output
phaseFFT = np.angle(calculatedFFT)      #calculates FFT phase from 
                                        #complex calculatedFFT output
phaseDegreesFFT = np.rad2deg(phaseFFT)  #convert to degrees

amplitudeScaledFFT = amplitudeFFT/float(dataLenghtFFT)
                 # scale by the number of points so that
                 # the magnitude does not depend on the length 
                 # of the signal
amplitudeScaledRMSFFT = amplitudeFFT/float(dataLenghtFFT)/math.sqrt(2)

# Scaling to Root mean square amplitude (dataLenghtFFT/sqrt{2}),
##############################################################################
# Plot the results
##############################################################################
plt.figure("FFT amplitude and phase coefficients")
plt.subplot(2,1,1)
plt.vlines(xFFT,0,amplitudeScaledFFT)
plt.title("FFT amplitude coefficients")
plt.xlabel("Harmonics")
plt.ylabel("Amplitude [V]")
plt.xlim(0,numberOfPoints/2+1) #adjuts the x axis to maximum of numberOfPoints
plt.grid(True)

plt.subplot(2,1,2)
plt.vlines(xFFT,0,phaseDegreesFFT)
plt.title("FFT phase coefficients")
plt.xlabel("Harmonics")
plt.ylabel("Phase [deg]")
plt.tight_layout()      #removes the overlapping of the labels in subplots
plt.xlim(0,numberOfPoints/2+1)
plt.grid(True)

plt.figure("Signal with fit and found minimum")
line, = plt.plot(np.rad2deg(x), y, linewidth=2)                 #fit
line, = plt.plot(np.rad2deg(x), ynoise, 'ro', linewidth=3)      #signal
plt.title("Signal with fit and found minimum")
plt.xlabel("Angle [deg]")
plt.ylabel("Amplitude [V]")
plt.xlim(0,numberOfPoints)
plt.vlines(FindMinimumPhase(phaseDegreesFFT[1]), 
           -amplitudeScaledFFT[1], amplitudeScaledFFT[1],linewidth=3,
            color='k',linestyles='dashed') 
            #draws the vertical line with the found minimum
plt.grid(True)
plt.show()
