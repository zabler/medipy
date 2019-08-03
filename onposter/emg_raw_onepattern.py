from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
EMG RAW ONE PATTERN

(1) Einlesen aller von PLUX26R in M6 und trennen nach EMG Signalen

(2) Signalwerte mit LsbValue und Baseline umrechnen

(3) Bereiche um Seizure 74910858 (PLUX26R) ausschneiden

(4) Signalstücke jeweils einzeln Plotten

Wiederholen für PLUX27L in M6 und Seizure 74911806

'''

# Input: Anzahl der Messungen, Signalart
anzahl_measurements = 2
channelliste1 = ['EMG1','EMG2','EMG3','EMG4','EMG7','EMG8']
labelliste1 = ['Pectoralis Major l','Neck extensors','Trizeps l','Biceps l','Deltoidus p. scapularis r', 'Trapezius p. descendenz r']
channelliste2 = ['EMG1','EMG2','EMG4','EMG5','EMG6','EMG7','EMG8']
labelliste2=['Sternocleidomastoideus r', 'Rectus Abdominis r', 'Rectus Femoris r', 'Pectoralis major r', 'Erector spinae','Biceps r', 'Triceps r']

# Fig für alle Signale
fig = plt.figure(frameon=True,figsize=(5,4))

# For Schleife über Anzahl der Messungen
for k in range(0, anzahl_measurements):
    if k == 0:
        # Movisensobjekt einlesen mit allen Signalarten einlesen
        movisensobject = m2pclass.m2pconverter(showtree=True)
        seizures = movisensobject.getentry('m6emgseizures').event
        # Nur die EMG Kanäle von PLUX26R nehmen
        for index, channelname in enumerate(channelliste1):
            # Signalparameter wählen und Signalwerte berechnen
            channel = movisensobject.getentry(channelname)
            channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
            fs = channel.sampleRate
            plt.plot(channel.signal[seizures[0]-20:seizures[0]+80], label=labelliste1[index],linewidth=0.7)
            
           
    elif k==1:
        # Movisensobjekt einlesen mit allen Signalarten einlesen
        movisensobject = m2pclass.m2pconverter(showtree=True)
        seizures = movisensobject.getentry('m6emgseizures').event
        # Nur die EMG Kanäle von PLUX227L nehmen
        for index, channelname in enumerate(channelliste2):
            # Signalparameter wählen und Signalwerte berechnen
            channel = movisensobject.getentry(channelname)
            channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
            fs = channel.sampleRate
            plt.plot(channel.signal[seizures[0]-20:seizures[0]+80], label=labelliste2[index],linewidth=0.7)

# Plot Seizure Onset            
plt.plot(20, 0, 'r--', label='seizure')
plt.axvline(x=20,color='r',linestyle='--')

# Plot Settings
#plt.title('Title',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('time [ms]',fontname="Arial")
plt.xlim(0, 100)
plt.ylim(-1.5,1.5)
plt.ylabel('EMG (different muscles) [mV]',fontname="Arial")
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=4)

# Beschriftung X-Achse neu
newtime = ['-20','0','20','40','60','80']
plt.gca().set_xticklabels(newtime)

# Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/emg_raw_onepattern.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/emg_raw_onepattern.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    
plt.show()

