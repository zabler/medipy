from wda.movisens2python import m2pclass
from wda.rdetections import dectclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
lookforHR
'''

# Input: Signalart    
channelname = 'ECG6'

# Movisensobjekt einlesen mit allen Signalarten einlesen
movisensobject = m2pclass.m2pconverter(showtree=True)

# Signalparameter wählen und Signalwerte berechnen
channel = movisensobject.getentry(channelname)
channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
fs = int(channel.sampleRate)

# Anfallsmarker einlesen
seizures = movisensobject.getentry('seizures').event
eR = int(movisensobject.getentry('seizures').sampleRate) #eventRate

# Frequenzverhältnis bestimmen
ver = fs/eR

#Abstand zum 
size=5*eR

#Timefaktor
timefaktor=60/((size/eR)*2)

#Vektor mit HR Berechnungszeitpunkten
#hrpoint = [-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60]
hrpoint=[]
for k in range(-21,20):
    hrpoint.append(k+1) 


for index, anfall in enumerate(seizures):
    hrwerte = []
    fig = plt.figure(index)
    for zeitpunkt in hrpoint:
            ana = channel.signal[int(((anfall+zeitpunkt*eR)-size)*ver):int(((anfall+zeitpunkt*eR)+size)*ver)]
            rpeaks_x, rr = dectclass.hamilton(ana, fs).detect()
            #rpeaks_x, rr = dectclass.pantompkins(ana, fs).detect()
            #rpeaks_x, rr = dectclass.ownr(ana, fs).detect()
            #rpeaks_x, rr = dectclass.skipi(ana, fs).detect()
            hrwerte.append((len(rpeaks_x)*timefaktor)) 
    #hrwerte = (hrwerte-(np.min(hrwerte)))/(np.max(hrwerte)-np.min(hrwerte))
    plt.plot(hrpoint,hrwerte,linewidth=0.5,color='#808080')
    plt.plot(hrpoint,hrwerte,".",color='#808080',label=f'Seizure #{index} Index {seizures[index]}')
    plt.axvline(x=0,color='r',linestyle='--')
    plt.plot(0, 0, 'r--', label='seizure onset')
    plt.ylim(min(hrwerte)-5,max(hrwerte)+5)
    plt.legend()
    
plt.show()
