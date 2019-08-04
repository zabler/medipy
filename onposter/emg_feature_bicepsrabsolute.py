from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import math


'''
EMG Feature Bicepsr Absolute

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

#Statistik
summe = np.zeros(100)
varliste =[]
counter=0

# For Schleife über Anzahl der Messungen

for k in range(0, anzahl_measurements):
    # Movisensobjekt einlesen mit allen Signalarten einlesen
    movisensobject = m2pclass.m2pconverter(showtree=True)
    #Anfallsmarker einlesen
    seizures = movisensobject.getentry('bicepsseizures').event
    eR = int(movisensobject.getentry('bicepsseizures').sampleRate)
    # Nur EMG8 des PLUX L nehmen
    channel = movisensobject.getentry(channelname)
    # Signalparameter wählen und Signalwerte berechnen
    fs = int(channel.sampleRate)
    channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
    #Frequenzverhätnis Samplerate zu Eventrate
    ver = fs/eR
    # Bereiche um Seizures ausschneiden
    size = 0.02 * eR #Size in Sekunden
    # For Schleife durch Anfälle
    for anfall in seizures:
        plotsignal = channel.signal[int((anfall-size)*ver):int((anfall+4*size)*ver)]
        #plotsignal = (plotsignal-(np.min(plotsignal)))/(np.max(plotsignal)-np.min(plotsignal)) #Comment for Absolut Average
        varliste.append(list(plotsignal))
        summe = summe + plotsignal
        #plt.plot(plotsignal,linewidth=0.5,color='#808080') #label=labelliste1[index])
        counter+=1

# Mittelwertbildung
mittel = summe/counter

# Streuungsbildung
var = np.zeros(100)

for ind in range(0, len(var)):
        varsumme = 0
        for k in varliste:
                a = (k[ind]-mittel[ind])**2
                varsumme = varsumme + a 
        var[ind] = math.sqrt(varsumme/(counter-1)) # Bessels Correction

# Plot Mittelwert
plt.plot(mittel,color='black',linewidth=2.0,label=f'relative arithmetic mean')

# Plot Var
plt.plot(mittel+var,color='black',linestyle='--',linewidth=1.5,label=f'standard deviation')
plt.plot(mittel-var,color='black',linestyle='--',linewidth=1.5)

# Plot Label for eech single Signal 
#plt.plot(0,0,color='#808080',label=f'EMG of n={counter} myoclonic seizures')

# Plot for Seizure Onset
plt.plot(size*ver, 0, 'r--', label='seizure onset')
plt.axvline(x=size*ver,color='r',linestyle='--')

#Plot Settings
#plt.title('Biceps r on different seizrues',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('time [ms]',fontname="Arial")
plt.ylabel('EMG (biceps r) [mV]',fontname="Arial") # Change Name for Absolut
plt.xlim(0, 100)
#plt.ylim(0,1)
#plt.ylim(-0.5,0.3) # Comment for Absolut
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=2)

# Beschriftung X-Achse neu
newtime = ['-20','0','20','40','60','80']
plt.gca().set_xticklabels(newtime)
# newwhy = ['0','20','40','60','80','100']
# plt.gca().set_yticklabels(newwhy)

# Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/emg_feature_bicepsrabsolute.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/emg_feature_bicepsrabsolute.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    
plt.show()
