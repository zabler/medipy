from wda.movisens2python import m2pclass
from wda.rdetections import dectclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import math

'''
EEG FEATURE FFT

(0) Konfiguration Filter 40 Hz

(1) Einlesen EEG Signals

(2) Signalwerte mit LsbValue und Baseline umrechnen

(3) Bereiche jeweils +- 5 Sekunden um Seizure ausschneiden 

(4) HR berechnen

(5) HR plotten

Wiederholen für PLUX27L in M6 und Seizure 74911806

'''

# Filter Configuration.
fS = 1000  # Sampling rate.
fL = 45  # Cutoff frequency.
N = 461 # Filter length, must be odd.
# Compute sinc filter.
tp = np.sinc(2 * fL / fS * (np.arange(N) - (N - 1) / 2))
# Apply window.
tp *= np.blackman(N)
# Normalize to get unity gain.
tp /= np.sum(tp)


# Input: Anzahl der Messungen, Signalart
anzahl_measurements = 7 #7
channelname='EEG5'
seizureseeg = [] #nur gute eeg

# Figur für alle Signale
fig = plt.figure(frameon=True,figsize=(5,4))

# Statistik
summe = np.zeros(400)
varliste=[]
counter = 0

# For Schleife über Anzahl der Messungen

for k in range(0, anzahl_measurements):

    # Movisensobjekt einlesen mit allen Signalarten einlesen
    movisensobject = m2pclass.m2pconverter(showtree=True)

    # Signalparameter wählen und Signalwerte berechnen
    channel = movisensobject.getentry(channelname)
    channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
    channel.signal = np.convolve(channel.signal,tp)
    fs = int(channel.sampleRate)

    # Anfallsmarker einlesen
    seizures = movisensobject.getentry('seizures').event
    # Nur gute EEG
    # for a in seizureseeg:
    #     if a in seizures:
    #         continue
    #     else:
    #         seizures.remove(a)

    eR = int(movisensobject.getentry('seizures').sampleRate) #eventRate
    
    # Frequenzverhältnis bestimmen
    ver = fs/eR

    #Abstand zum Anfall in Sekunden
    size=5*eR


    # Calculation
    for anfall in seizures:
        ana = channel.signal[int((anfall-size)*ver):int((anfall+size)*ver)]
        freqwerte = np.fft.fft(ana)
        halbe= int(len(freqwerte)/2 +1)
        freqwerte = abs(2*freqwerte[:halbe]/halbe)
        frequenzen = np.linspace(0,fs/2,halbe)
        freqwerte = freqwerte[0:400]
        frequenzen = frequenzen[0:400]
        
        #freqwerte = (freqwerte-(np.min(freqwerte)))/(np.max(freqwerte)-np.min(freqwerte))
        varliste.append(freqwerte)
        summe += freqwerte
        counter+=1
        plt.plot(frequenzen,freqwerte,linewidth=0.5,color='#808080')
  

# Mittelwertbildung
mittel = summe/counter

# Streuungsbildung
var = np.zeros(400)

for ind in range(0, len(var)):
        varsumme = 0
        for k in varliste:
                a = (k[ind]-mittel[ind])**2
                varsumme = varsumme + a 
        var[ind] = math.sqrt(varsumme/(counter-1)) # Bessels Correction

#Plot Mittelwert
plt.plot(frequenzen,mittel,color='black',linewidth=2.0,label=f'arithmetic mean')
#plt.plot(hrpoint,mittel,".",color='black')

# Plot Var
plt.plot(frequenzen,mittel+var,color='black',linestyle='--',linewidth=1.5,label=f'standard deviation')
plt.plot(frequenzen,mittel-var,color='black',linestyle='--',linewidth=1.5)

#Plot for Labeling single Signals
plt.plot(0,0,color='#808080',label=f'EEG FFTs of n={counter} myoclonic seizures')


#Plot Settings
#plt.title('Biceps r on different seizrues',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('Frequency [Hz]',fontname="Arial")
plt.ylabel('Magnitude',fontname="Arial")
plt.xlim(0, 40)
#plt.ylim(0,30)
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=2)

# Beschriftung X-Achse neu
# newtime = ['-180','-60','','','','10','30','60','180']
# plt.gca().set_xticks([-180,-60,-30,-10,0,10,30,60,180])
# plt.gca().set_xticklabels(newtime)
# newwhy = ['0','20','40','60','80','100']
# plt.gca().set_yticklabels(newwhy)


#Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/eeg_feature_freq.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/eeg_feature_freq.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    

plt.show()
