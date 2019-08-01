from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


'''
EMG ALL BICEPS ONE

(1) Einlesen EMG Signal Biceps r und Anfallsmarker

(2) Signalwerte mit LsbValue umrechnen

(3) Bereiche der Signale um Seizures nehmen, auf % Spreizen und zur Summe hinzuaddieren, Plotten

(4) Wiederholen für alle Messreihen

(5) Mittelwert bilden und Plotten

'''

# Input: Anzahl der Messungen, Signalart
anzahl_measurements = 7
channelname='EMG8'
labelliste1 = ['Biceps r']


# Fig für alle Signale
fig = plt.figure(frameon=True,figsize=(5,4))

#Mittelwertbilder
summe = np.zeros(100)
counter=0

# For Schleife über Anzahl der Messungen

for k in range(0, anzahl_measurements):
    # Movisensobjekt einlesen mit allen Signalarten einlesen
    movisensobject = m2pclass.m2pconverter(showtree=True)
    #Anfallsmarker einlesen
    seizures = movisensobject.getentry('bicepsseizures').event
    # Nur EMG8 des PLUX L nehmen
    channel = movisensobject.getentry(channelname)
    # Signalparameter wählen und Signalwerte berechnen
    fs = channel.sampleRate
    channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
    # For Schleife durch Anfälle
    for anfall in seizures:
        plotsignal = channel.signal[anfall-20:anfall+80]
        plotsignal = (plotsignal-(np.min(plotsignal)))/(np.max(plotsignal)-np.min(plotsignal)) #Comment for Absolut Average
        summe = summe + plotsignal
        plt.plot(plotsignal,linewidth=0.5,color='#808080') #label=labelliste1[index])
        counter+=1

# Mittelwertbildung
summe = summe/counter

# Plot Mittelwert
plt.plot(summe,color='black',linewidth=2.0,label=f'relative average emg ictal')

# Plot Label for eech single Signal 
plt.plot(0,0,color='#808080',label=f'spread emg voltage of n={counter} myoclonic seizures')

# Plot for Seizure Onset
plt.plot(20, 0, 'r--', label='Onset')
plt.axvline(x=20,color='r',linestyle='--')

#Plot Settings
#plt.title('Biceps r on different seizrues',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('time [ms]',fontname="Arial")
plt.ylabel('equalized emg voltage (biceps r) [%]',fontname="Arial") # Change Name for Absolut
plt.xlim(0, 100)
#plt.ylim(-1.5,1.5) # Comment for Absolut
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=2)


# Beschriftung X-Achse neu
newtime = ['-20','0','20','40','60','80']
plt.gca().set_xticklabels(newtime)

# Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/emg_biceps.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/emg_biceps.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    
plt.show()
