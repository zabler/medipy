from __future__ import division
from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import math

'''
EDA FEATURE ABSOLUT

(0) Defintion Tiefpassfilter

(1) Einlesen eines eda Signals und Anfallsmarker

(2) Signalwerte mit LsbValue umrechnen

(3) Bereiche der Signale um Seizures nehmen und zur Summe hinzuaddieren

(4) Wiederholen für alle Messreihen

(5) Mittelwert bilden und Plotten

'''

# Filter Configuration.
fS = 32  # Sampling rate.
fL = 2  # Cutoff frequency.
N = 47 # Filter length, must be odd.
# Compute sinc filter.
tp = np.sinc(2 * fL / fS * (np.arange(N) - (N - 1) / 2))
# Apply window.
tp *= np.blackman(N)
# Normalize to get unity gain.
tp /= np.sum(tp)

# Input: Anzahl der Messungen, Signalart
anzahl_measurements = 7
channelname='eda'

# Fig für alle Signale
fig = plt.figure(frameon=True,figsize=(5,4))

# Statistik
summe = np.zeros(800)
varliste=[]
counter = 0

# For Schleife über Anzahl der Messungen
for k in range(0, anzahl_measurements):
    # Movisensobjekt einlesen mit allen Signalarten einlesen
    movisensobject = m2pclass.m2pconverter(showtree=True)
    # Signalparameter wählen und Signalwerte berechnen
    channel = movisensobject.getentry(channelname)
    channel.signal = channel.signal * float(channel.lsbValue)
    fs = int(channel.sampleRate)
    # Tiefpassfiltern
    channel.signal = np.convolve(channel.signal,tp)
    #Anfallsmarker einlesen
    seizures = movisensobject.getentry('seizures').event
    eR = int(movisensobject.getentry('seizures').sampleRate)
    # Frequenzverhältnis
    ver = fs/eR
    #Abstand zum 
    size=5*eR
    for anfall in seizures:
        plotsignal = channel.signal[int((anfall-size)*ver):int((anfall+4*size)*ver)]
        varliste.append(list(plotsignal))
        summe += plotsignal
        counter +=1
        #plt.plot(plotsignal,linewidth=0.5,color='#808080') #Uncomment for Single Plots

#Mittelwertbildung
mittel = summe/counter

# Streuungsbildung
var = np.zeros(800)

for ind in range(0, len(var)):
        varsumme = 0
        for k in varliste:
                a = (k[ind]-mittel[ind])**2
                varsumme = varsumme + a 
        var[ind] = math.sqrt(varsumme/(counter-1)) # Bessels Correction

#Plot Mittelwert
plt.plot(mittel,color='black',linewidth=2.0,label=f'relative arithmetic mean')

# Plot Var
plt.plot(mittel+var,color='black',linestyle='--',linewidth=1.5,label=f'standard deviation')
plt.plot(mittel-var,color='black',linestyle='--',linewidth=1.5)


# #Plot for Labeling single Signals #Uncoment for Single Plots
# plt.plot(0,0,color='#808080',label=f'absoult SC of n={counter} myoclonic seizures')

# Seizure Onset plotten
plt.plot(size*ver, 0, 'r--', label='seizure onset')
plt.axvline(x=size*ver,color='r',linestyle='--')

# Settings Plot
#plt.title('Title',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('time [s]',fontname="Arial")
plt.xlim(0, 800)
#plt.ylim(4,7)
plt.ylabel('skin conductance (palm) [µS]',fontname="Arial")
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=2)

# Beschriftung X-Achse neu
newtime = ['-5','0','5','10','15','20']
plt.gca().set_xticks([0,160,320,480,640,800])
plt.gca().set_xticklabels(newtime)

# Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/eda_feature_absolut.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/eda_feature_absolut.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    

plt.show()
