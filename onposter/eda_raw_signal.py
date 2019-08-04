from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
EDA RAW SIGNAL

(1) Erstellung eines Movisens Objekts mit ACM Signal

(2) Signalwerte mit LsbValue umrechnen

(3) Bereich um Seizure berechnen

(4) Signalstück einzeln Plotten

'''

# Input: Signalart
channelname = 'eda'

# Figur erstellen
fig = plt.figure()

# Movisensobjekt einlesen mit gewünschter Signalart 
movisensobject = m2pclass.m2pconverter(
    channelname, 'seizures', showtree=True)

# SignalEntry ECG und ValuesEntry NN_Live auslesen
channel = movisensobject.getentry(channelname)
fs = int(channel.sampleRate)
channel.signal = channel.signal * float(channel.lsbValue)

# Seizures: Anfall wählen
seizures = movisensobject.getentry('seizures').event
anfall = seizures[4] #Seizure 6 in M3
eR = int(movisensobject.getentry('seizures').sampleRate)

#Frequenzverhätnis Samplerate zu Eventrate
ver = fs/eR
size = 5 * eR #Size in Sekunden

# Bereiche um Seizures ausschneiden
# x samples im Signal mit fs Abtastfrequenz = (anfall-size) *ver = (anfall*eR - T*eR) * fs/eR = (anfall-T)*fs
plt.plot(channel.signal[int((anfall-size)*ver):int((anfall+4*size)*ver)], label='skin conductance', linewidth=0.7, color='black')

# Seizure Onset plotten
plt.plot(size*ver, 0, 'r--', label='seizure onset')
plt.axvline(x=size*ver,color='r',linestyle='--')

# Settings Plot
#plt.title('Title',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('time [s]',fontname="Arial")
plt.xlim(0, 800)
plt.ylim(2,5)
plt.ylabel('skin conductance (thenar) [µS]',fontname="Arial")
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=2)

# Beschriftung X-Achse neu
newtime = ['-5','0','5','10','15','20']
plt.gca().set_xticks([0,160,320,480,640,800])
plt.gca().set_xticklabels(newtime)

# Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/eda_raw_signal.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/eda_raw_signal.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    

plt.show()
