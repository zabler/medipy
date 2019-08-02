from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
EEG Raw t5p3

(1) Erstellung eines Movisens Objekts mit ACM Signal

(2) Signalwerte mit LsbValue umrechnen

(3) Bereich um Seizure berechnen
M1 EEG5 #4 
M1 EEG5 #8

(4) Signalstück einzeln Plotten

'''

# Input: Signalart
channelname = 'EEG5'

# Figur erstellen
fig = plt.figure()

# Movisensobjekt einlesen mit gewünschter Signalart 
movisensobject = m2pclass.m2pconverter(
    channelname, 'seizures', showtree=True)

# SignalEntry ECG und ValuesEntry NN_Live auslesen
channel = movisensobject.getentry(channelname)
fs = int(channel.sampleRate)
channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)

# Seizures: Anfall wählen
seizures = movisensobject.getentry('seizures').event
anfall = seizures[4] # Seizure 4 in M1
eR = int(movisensobject.getentry('seizures').sampleRate)

#Frequenzverhätnis Samplerate zu Eventrate
ver = fs/eR

# Bereiche um Seizures ausschneiden
size = 1 * eR #Size in Sekunden

# Bereiche um Seizures ausschneiden
# x samples im Signal mit fs Abtastfrequenz = (anfall-size) *ver = (anfall*eR - T*eR) * fs/eR = (anfall-T)*fs
plt.plot(channel.signal[int((anfall-size)*ver):int((anfall+2*size)*ver)], label='EEG T5-P3 ictal', linewidth=0.7, color='black')
   
# Plot Seizure Onset
plt.plot(size*ver, 0, 'r--', label='seizure onset')
plt.axvline(x=size*ver,color='r',linestyle='--')

# Plot Settings
#plt.title('Ttile',fontname="Arial", fontweight="bold",loc='left')
plt.xlabel('time [ms]',fontname="Arial")
plt.xlim(0, 3000)
plt.ylim(-40,40)
plt.ylabel('EEG (T5-P3) [µV]',fontname="Arial")
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=4)

#Beschriftung X-Achse neu
newtime = ['-1000','-500','0','500','1000','1500','2000']
#plt.gca().set_xticks([0,16,32,48,64])
plt.gca().set_xticklabels(newtime)

# Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/eeg_raw_t5p3.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/eeg_raw_t5p3.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    

plt.show()
