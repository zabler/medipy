from wda.movisens2python import m2pclass
from wda.rdetections import dectclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
ECG Raw Artifact

(1) Erstellung eines Movisens Objekts mit ACM Signal

(2) Signalwerte mit LsbValue umrechnen

(3) Bereich um Seizure berechnen

(4) Signalstück einzeln Plotten

'''

# Input: Signalart
channelname = 'ecg'

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
anfall = seizures[7] #Seizure 7 in M6
eR = int(movisensobject.getentry('seizures').sampleRate)

#Frequenzverhätnis Samplerate zu Eventrate
ver = fs/eR
size = 0.5 * eR #Size in Sekunden

# Bereiche um Seizures ausschneiden
# x samples im Signal mit fs Abtastfrequenz = (anfall-size) *ver = (anfall*eR - T*eR) * fs/eR = (anfall-T)*fs
ana = channel.signal[int((anfall-size)*ver):int((anfall+2*size)*ver)]
rpeaks_x, rr = dectclass.skipi(ana, fs).detect()
plt.plot(channel.signal[int((anfall-size)*ver):int((anfall+2*size)*ver)], label='ECG with seizure artifact', linewidth=0.7, color='black')
plt.plot(rpeaks_x,ana[rpeaks_x],'x',color='red',label='detected R-Peaks')   

# Plot Seizure Onset
plt.plot(size*ver, 0, 'r--', label='seizure onset')
plt.axvline(x=size*ver,color='r',linestyle='--')

# Plot Settings
#plt.title('Title',fontname="Arial", fontweight="bold",loc='left')
plt.xlabel('time [ms]',fontname="Arial")
plt.xlim(0, 384)
plt.ylim(-1,3)
plt.ylabel('ECG (chest) [mV]',fontname="Arial")
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=4)

# Beschriftung X-Achse neu
newtime = ['-500','-250','0','250','500','750','1000']
plt.gca().set_xticks([0,64,128,192,256,320,384])
plt.gca().set_xticklabels(newtime)

# Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/ecg_raw_artifact.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/ecg_raw_artifact.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    

plt.show()
