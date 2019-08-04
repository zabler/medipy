from wda.movisens2python import m2pclass
from wda.rdetections import dectclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import math

'''
ECG FEATURE HR ASOLUTE

(1) Einlesen ECG Signal

(2) Signalwerte mit LsbValue und Baseline umrechnen

(3) Bereiche jeweils +- 15 Sekunden um Seizure ausschneiden 
-1min, -30, -10, 0, +10, +30, +1min, +5min

(4) HR berechnen

(5) HR plotten

Wiederholen für PLUX27L in M6 und Seizure 74911806

'''

# Input: Anzahl der Messungen, Signalart
anzahl_measurements = 7
channelname='ecg'

# Figur für alle Signale
fig = plt.figure(frameon=True,figsize=(5,4))

# Statistik
summe = np.zeros(121)
varliste=[]
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
    
    # Frequenzverhältnis bestimmen
    ver = fs/eR

    #Abstand zum 
    size=5*eR

    #Vektor mit HR Berechnungszeitpunkten
    hrpoint=[]
    for k in range(-61,60):
         hrpoint.append(k+1) 
    
    #Timefaktor
    timefaktor = 60/(5*2) #Change Size faktor


    # Calculation
    for anfall in seizures:
        hrwerte = []
        for zeitpunkt in hrpoint:
                ana = channel.signal[int(((anfall+zeitpunkt*eR)-size)*ver):int(((anfall+zeitpunkt*eR)+size)*ver)]
                #rpeaks_x, rr = dectclass.hamilton(ana, fs).detect()
                rpeaks_x, rr = dectclass.pantompkins(ana, fs).detect()
                #rpeaks_x, rr = dectclass.ownr(ana, fs).detect()
                #rpeaks_x, rr = dectclass.skipi(ana, fs).detect()
                hrwerte.append((len(rpeaks_x)*timefaktor)) #BPM = Anzahl Rzacken geteilt durch 2*Sizefaktor mal 60  
        #hrwerte = (hrwerte-(np.min(hrwerte)))/(np.max(hrwerte)-np.min(hrwerte))
        varliste.append(hrwerte)
        summe += hrwerte
        counter+=1
        #plt.plot(hrpoint,hrwerte,linewidth=0.5,color='#808080')
        #plt.plot(hrpoint,hrwerte,".",color='#808080')
       

# Mittelwertbildung
mittel = summe/counter

# Streuungsbildung
var = np.zeros(121)

for ind in range(0, len(var)):
        varsumme = 0
        for k in varliste:
                a = (k[ind]-mittel[ind])**2
                varsumme = varsumme + a 
        var[ind] = math.sqrt(varsumme/(counter-1)) # Bessels Correction

#Plot Mittelwert
plt.plot(hrpoint,mittel,color='black',linewidth=2.0,label=f'arithmetic mean')
#plt.plot(hrpoint,mittel,".",color='black')

# Plot Var
plt.plot(hrpoint,mittel+var,color='black',linestyle='--',linewidth=1.5,label=f'standard deviation')
plt.plot(hrpoint,mittel-var,color='black',linestyle='--',linewidth=1.5)

#Plot for Labeling single Signals
#plt.plot(0,0,color='#808080',label=f'HR of n={counter} myoclonic seizures')

# #Plot Seizure Onset
plt.plot(0, 0, 'r--', label='seizure onset')
plt.axvline(x=0,color='r',linestyle='--') #Osec bei Onset

#Plot Settings
#plt.title('Biceps r on different seizrues',fontname="Arial", fontweight="bold",loc='left') #fontsize=12
plt.xlabel('time [s]',fontname="Arial")
plt.ylabel('HR [bpm]',fontname="Arial")
plt.xlim(-60, 60)
#plt.ylim(40,150) with single plots
plt.ylim(60,120)
plt.grid(b=True,which='major',axis='both')
plt.legend(fontsize='xx-small',bbox_to_anchor=(0,1.02,1,0.5), loc="lower left",mode='expand',borderaxespad=0, ncol=2)

# Beschriftung X-Achse neu
# newtime = ['-180','-60','','','','10','30','60','180']
# plt.gca().set_xticks([-180,-60,-30,-10,0,10,30,60,180])
# plt.gca().set_xticklabels(newtime)
# newwhy = ['0','20','40','60','80','100']
# plt.gca().set_yticklabels(newwhy)


#Bilder speichern
plt.savefig('/Users/nicolaszabler/Desktop/ecg_feature_hrabsolute.png',dpi=300,transparent=False,bbox_inches='tight')    
plt.savefig('/Users/nicolaszabler/Desktop/ecg_feature_hrabsolute.svg',dpi=300,format='svg',transparent=False, bbox_inches='tight')    

plt.show()
