from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
# from matplotlib import style
# style.use('fivethirtyeight')
import matplotlib
matplotlib.use('TkAgg')

'''
EMG ALL 

(1) Einlesen aller von PLUX26R in M6 und trennen nach EMG Signalen

(2) Signalwerte mit LsbValue und Baseline umrechnen

(3) Bereiche um Seizure 74910858 (PLUX26R) ausschneiden

(4) Signalstücke jeweils einzeln Plotten

Wiederholen für PLUX27L in M6 und Seizure 74911806

'''

# Input: Anzahl der Messungen, Signalart
anzahl_measurements = 6
channelname='EMG8'
labelliste1 = ['Biceps r']


# Fig für alle Signale
fig = plt.figure(frameon=True,figsize=(5,4))
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
    #for index, anfall in enumerate(seizures):
    for anfall in seizures:
        plotsignal = channel.signal[anfall-20:anfall+80]
        plotsignal = (plotsignal-(np.min(plotsignal)))/(np.max(plotsignal)-np.min(plotsignal))
        summe = summe + plotsignal
        plt.plot(plotsignal,linewidth=0.5,color='#808080') #label=labelliste1[index])
        counter+=1

summe = summe/counter
labelname = f'relative average emg ictal'
labelname2 = f'spread emg voltage of n={counter} myoclonic seizures'
plt.plot(summe,color='black',linewidth=2.0,label=labelname)
plt.plot(0,0,color='#808080',label=labelname2)
plt.plot(20, 0, 'r--', label='Onset')
plt.axvline(x=20,color='r',linestyle='--')
#plt.title('Biceps r on different seizrues',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('time [ms]',fontname="Arial")
plt.xlim(0, 100)
plt.ylim(-1.5,1.5)
#plt.ylabel('emg voltage (biceps r) [mV]',fontname="Arial")
plt.ylabel('equalized emg voltage (biceps r) [%]',fontname="Arial")
plt.grid(b=True,which='major',axis='both')
#plt.legend(fontsize='xx-small',bbox_to_anchor=(0,-0.52,1,0.4), loc="upper left",mode='expand',borderaxespad=0, ncol=7)
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=4)
newtime = ['-20','0','20','40','60','80']
plt.gca().set_xticklabels(newtime)
plt.savefig('/Users/nicolaszabler/Desktop/emg_biceps.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/emg_biceps.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    
plt.show()
