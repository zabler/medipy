from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
ACHq<ey<s

(1) Erstellung eines Movisens Objekts mit ACM Signal

(2) Signalwerte mit LsbValue umrechnen

(3) Bereich um Seizure berechnen

(4) Signalstücke der Subchannels jeweils einzeln Plotten

'''

# Input: Signalart
channelname = 'acc'

# Figur erstellen
fig = plt.figure()

# For Schleife über Anzahl der Messungen

# Movisensobjekt einlesen mit gewünschter Signalart 
movisensobject = m2pclass.m2pconverter(
    channelname, 'seizures', showtree=True)

# Signalparameter wählen und Signalwerte berechnen
channel = movisensobject.getentry(channelname)
fs = channel.sampleRate
le = len(channel.signal)

# Subchannels auslesen
accX = channel.signal[0:le:3] * float(channel.lsbValue)
accY = channel.signal[1:le:3] * float(channel.lsbValue)
accZ = channel.signal[2:le:3] * float(channel.lsbValue)

# Marker wählen / Freqverhältnis Samplerate und Eventrate ist bei beiden gleich
seizures = movisensobject.getentry('seizures').event
anfall = seizures[4] #Seizure 3533865 in M3

# Bereiche um Seizures ausschneiden
plt.plot(accX[anfall-16:anfall+48], label='ACC X', linewidth=0.7, color='black')
plt.plot(accY[anfall-16:anfall+48], label='ACC Y', linewidth=0.7, color='blue')
plt.plot(accZ[anfall-16:anfall+48], label='ACC Z', linewidth=0.7, color='green')
   

# Plot Seizure Onset
plt.plot(16, 0, 'r--', label='Onset')
plt.axvline(x=16,color='r',linestyle='--')

# Plot Settings
#plt.title('Different EMGs of one Seizure',fontname="Arial", fontweight="bold",loc='left')
plt.xlabel('time [ms]',fontname="Arial")
plt.xlim(0, 64)
plt.ylim(0,1)
plt.ylabel('3D ACC (chest) [g]',fontname="Arial")
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=4)

# Beschriftung X-Achse neu
newtime = ['-250','0','250','500','750']
plt.gca().set_xticks([0,16,32,48,64])
plt.gca().set_xticklabels(newtime)

# Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/acm_chest.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/acm_chest.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    

plt.show()
