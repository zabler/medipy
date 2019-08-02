from wda.movisens2python import m2pclass
from wda.rdetections import dectclass
import numpy as np
from matplotlib import pyplot as plt
# from matplotlib import style
# style.use('fivethirtyeight')
import matplotlib
matplotlib.use('TkAgg')
import math

'''
EMG ALL 

(1) Einlesen aller von PLUX26R in M6 und trennen nach EMG Signalen

(2) Signalwerte mit LsbValue und Baseline umrechnen

(3) Bereiche um Seizure 74910858 (PLUX26R) ausschneiden

(4) Signalstücke jeweils einzeln Plotten

Wiederholen für PLUX27L in M6 und Seizure 74911806

'''

# Input: Anzahl der Messungen, Signalart
anzahl_measurements = 1 #6
channelname='ecg'

# Figur für alle Signale
fig = plt.figure(frameon=True,figsize=(5,4))
counter = 0

# For Schleife über Anzahl der Messungen

for k in range(0, anzahl_measurements):
    # Movisensobjekt einlesen mit allen Signalarten einlesen
    movisensobject = m2pclass.m2pconverter(showtree=True)
    # Signalparameter wählen und Signalwerte berechnen
    channel = movisensobject.getentry(channelname)
    channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
    fs = int(channel.sampleRate)
    # Anfallsmarker einlesen
    seizures = movisensobject.getentry('seizures').event
    eR = int(movisensobject.getentry('seizures').sampleRate) #eventRate
    ver = fs/eR
    for anfall in seizures:
        ana = channel.signal[math.ceil(anfall*ver)-30720:math.ceil(anfall*ver)+30720]
        rpeaks_x, rr = dectclass.hamilton(ana, fs).detect()
        counter +=1
        plt.plot(rr,linewidth=0.5,color='#808080')
       

#summe = summe/counter
labelname = f'average RR interval'
labelname2 = f'Ictal RR intervals of n={counter} myoclonic seizures'
#plt.plot(summe,color='black',linewidth=2.0,label=labelname)
# plt.plot(0,0,color='#808080',label=labelname2)
# plt.plot(120, 0, 'r--', label='Onset')
# plt.axvline(x=120,color='r',linestyle='--')
#plt.title('Biceps r on different seizrues',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('time [min]',fontname="Arial")
#plt.xlim(0, 240)
#plt.ylim(0,50)
plt.ylabel('RR interval [ms]',fontname="Arial")
plt.grid(b=True,which='major',axis='both')
#plt.legend(fontsize='xx-small',bbox_to_anchor=(0,-0.52,1,0.4), loc="upper left",mode='expand',borderaxespad=0, ncol=7)
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=4)
# newtime = ['-2','','','','0','','','','2']
# plt.gca().set_xticklabels(newtime)
# plt.savefig('/Users/nicolaszabler/Desktop/eda_all.png',dpi=300,transparent=False,bbox_inches='tight')    
# plt.savefig('/Users/nicolaszabler/Desktop/eda_all.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    
plt.show()




       
            # if k in timestamps:
            #     plotsignal.append(werte[timestamps.index(k)])
            # elif len(plotsignal) == 0:
            #     continue
            # else:
            #     plotsignal.append(plotsignal[-1])
           # Globale Variablen zur Mittelwertbildung
# summe = np.zeros(122880)
# counter = 0
 #plt.plot(channel.signal[math.ceil(anfall*ver)-160:math.ceil(anfall*ver)+640]) #ZEIT #porzenutal shiten?
        #plotsignal = (plotsignal-(np.min(plotsignal)))/(np.max(plotsignal)-np.min(plotsignal))
        #summe += plotsignal
