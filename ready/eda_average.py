from __future__ import division
from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import math

'''
EDA AVERAGE

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

# Mittelwertbilder
summe = np.zeros(800)
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
    # Frequenzverhältnis von Samplerate zu Eventrate
    ver = fs/eR
    for anfall in seizures:
        plotsignal = channel.signal[math.ceil(anfall*ver)-160:math.ceil(anfall*ver)+640]
        #plotsignal = (plotsignal-(np.min(plotsignal)))/(np.max(plotsignal)-np.min(plotsignal))
        summe += plotsignal
        counter +=1

# Mittelwertbildung
summe = summe/counter

# Mittelwert plotten
plt.plot(summe,color='black',linewidth=3.0,label=f'absolute average SC post-ictal (n={counter})')

# Seizure Onset plotten
plt.plot(160, 0, 'r--', label='Onset')
plt.axvline(x=160,color='r',linestyle='--')

# Settings Plot
#plt.title('Biceps r on different seizrues',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('time [s]',fontname="Arial")
plt.xlim(0, 800)
plt.ylim(4,7)
plt.ylabel('skin conductance (palm) [µS]',fontname="Arial")
plt.grid(b=True,which='major',axis='both')
#plt.legend(fontsize='xx-small',bbox_to_anchor=(0,-0.52,1,0.4), loc="upper left",mode='expand',borderaxespad=0, ncol=7)
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=4)

# Beschriftung X-Achse neu
newtime = ['-5','','','','12.5','','','','20']
plt.gca().set_xticklabels(newtime)

# Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/eda_average.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/eda_average.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    

plt.show()
