import csv
import struct
import os
import numpy as np
import pandas as pd
import xml.etree.ElementTree as et
from tkinter import filedialog, Tk
import time

# Only for Example
import numpy as np
from matplotlib import pyplot as plt

'''BASECLASS'''


class movisens2python():
    '''
    Klasse zum einlesen und bereitstellen von Movisensdaten in Python, zuvor abgespeichert im CSV-Format
    '''
    # Verbesserungsmöglichkeiten:
    # - Daten schön machen inkl. Datasheet infos z.B CM = 560µV ECG
    # - NoData Option: Signale sollen nicht eingelesen werden aber dafür XML Info,
    #   dazuSignal getrennt einlesen, nachdem XML Info geladen bereis fertig geladen wurde: Extra Funktionen für Read
    #   dazu nötig: Subclassfunktion die den Parameter bearbeitet
    # - Default Ausgabe: Baumstruktur nur auch vorhandene Parameter, bei den Klassen werden alle Attribute angezeigt auch wenn nicht vorhanden
    # - Allgemeins Problem: Mögliche Leere Auswahlparameter bei verschiedenen Entrys
    # - z.B Evententry hat kein Attr. Signal, und nicht jeder Signalentry hat eine Baseline
    # - Zeitbereichsauswahl ermöglichen: Start und Dauer
    # - @propertys einführen für choosedata und customsettings
    # - Exits und Fehlermeldungen einbauen

    def __init__(self):
        '''Init der BaseClass Objektvariablen'''
        # Choosedata
        self.filename = None
        self.filepfad = None
        self.name = None
        # Customsettings
        self.start = None
        self.dauer = None
        self.readinlist = []
        self.defaultcontents = ['acc', 'angularrate', 'artifact', 'bpmbxb_live', 'charging', 'ecg', 'hr_live', 'hrvisvalid_live', 'hrvrmssd_live', 'marker',
                                'movementacceleration_live', 'nn_live','nnlist', 'press', 'stateofcharge', 'stepcount_live', 'temp', 'tempmean_live','eda',
                                'EMG1','EMG2','EMG3','EMG4','EMG5','EMG6','EMG7','EMG8','EEG3','EEG5','ECG6','seizures','m6seizures','m6emgseizures','bicepsseizures']

    def movisens_choosedata(self, extension='.xml', datatype='unisens'):
        '''Funktion zum Auswählen der Daten'''
        root = Tk()  # TinkerDialog öffnen
        root.withdraw()  # Tinkerfenster verstecken
        files = filedialog.askopenfilenames(
            defaultextension='.*',
            filetypes=[(datatype, extension), ],
            title="Select Files",
            multiple=True,
        )
        self.filepfad = os.path.dirname(files[0])
        if len(files) <= 1:
            self.filename = os.path.split(files[0])[1]
        else:
            self.filename = os.path.split(self.filepfad)[1]
        self.name = os.path.splitext(self.filename)[0]

    def movisens_customsettings(self, Datenliste, Start, Dauer):
        '''Funktion zur Bearbeitung der Custom Input Parameter'''
        if not Datenliste:
            self.readinlist = self.defaultcontents
        else:
            self.readinlist = Datenliste

        self.start = Start
        self.dauer = Dauer

    def movisens_add_content(self, name):
        self.defaultcontents.append(name)

    def movisens_get_XML(self):
        '''Funktion zum Lesen der XML-Unisens Informationen '''

        # Namespace und Wurzel deklarieren
        unisensspace = {'uni': 'http://www.unisens.org/unisens2.0'}
        tree = et.parse(self.filepfad + '/' + self.filename)
        wurzel = tree.getroot()

        # RootAttributes
        self.rootinfo = rootAttributes()
        for key, val in wurzel.attrib.items():
            if hasattr(self.rootinfo, key):
                setattr(self.rootinfo, key, val)

        # CustomAttributes
        self.custominfo = customAttributes()
        for customatt in wurzel.findall('uni:customAttributes/uni:customAttribute', namespaces=unisensspace):
            if hasattr(self.custominfo, customatt.attrib['key']):
                setattr(
                    self.custominfo, customatt.attrib['key'], customatt.attrib['value'])

        # SignalEntry
        for signal in wurzel.findall('uni:signalEntry', namespaces=unisensspace):
            if signal.attrib['id'].split('.')[0] in self.readinlist:

                # Signaleintrag erstellen
                self.signalEntry = signalEntry()

                # Attributes speichern
                for key, value in signal.attrib.items():
                    if hasattr(self.signalEntry, key):
                        setattr(self.signalEntry, key, value)

                # Channels speichern
                channellist = []
                for channel in signal.findall('uni:channel', namespaces=unisensspace):
                    channellist.append(channel.attrib['name'])
                setattr(self.signalEntry, 'channels', channellist)

                # Signaldatei einlesen
                # Prüfen ob Binary oder CSV
                if signal.attrib['id'].split('.')[1] == 'csv':
                    with open(str(self.filepfad + '/' + signal.attrib['id']), 'r') as csvfile:
                        read = np.genfromtxt(
                            csvfile, dtype=np.int, delimiter=';')  # csv nimmt automatisch richtige integersize
                        # damit alle gleiches Format haben
                        read = read.astype('int32')
                else:  # hier nochmal unterscheiden zwischen int16 und int32
                    if signal.attrib['dataType'] == 'int32':
                        with open(str(self.filepfad + '/' + signal.attrib['id']), 'rb') as binfile:
                            read = np.fromfile(binfile, dtype=np.int32)
                    else:
                        with open(str(self.filepfad + '/' + signal.attrib['id']), 'rb') as binfile:
                            read = np.fromfile(binfile, dtype=np.int16)
                            # damit alle gleiches Format haben
                            read = read.astype('int32')

                setattr(self.signalEntry, 'signal', read)

                # SignalEntry in Namen der jeweiligen CSV Datei umschreiben
                rename_attribute(self, 'signalEntry',
                                 signal.attrib['id'].split('.')[0])

        # ValuesEntry
        for values in wurzel.findall('uni:valuesEntry', namespaces=unisensspace):
            if values.attrib['id'].split('.')[0] in self.readinlist:

                # Valueseintrag erstellen
                self.valuesEntry = valuesEntry()

                # Attributes speichern
                for key, value in values.attrib.items():
                    if hasattr(self.valuesEntry, key):
                        setattr(self.valuesEntry, key, value)

                # Channels speichern
                channellist = []
                for channel in values.findall('uni:channel', namespaces=unisensspace):
                    channellist.append(channel.attrib['name'])
                setattr(self.valuesEntry, 'channels', channellist)

                # Valuesdatei einlesen
                # If not needed Values Entry only as CSV
                if values.attrib['id'].split('.')[1] == 'csv':
                    with open(str(self.filepfad + '/' + values.attrib['id']), 'r') as csvfile:
                        read = np.genfromtxt(
                            csvfile, dtype=np.uint, delimiter=';')
                        # read = list(csv.reader(csvfile, delimiter=';'))
                        # read = np.array(read[0:], dtype=np.int)

                setattr(self.valuesEntry, 'values', read)

                # ValuesEntry in Namen der jeweiligen CSV Datei umschreiben
                rename_attribute(self, 'valuesEntry',
                                 values.attrib['id'].split('.')[0])

        # EventEntry
        for event in wurzel.findall('uni:eventEntry', namespaces=unisensspace):
            if event.attrib['id'].split('.')[0] in self.readinlist:

                # Eventeintrag erstellen
                self.eventEntry = eventEntry()

                # Attributes speichern
                for key, value in event.attrib.items():
                    if hasattr(self.eventEntry, key):
                        setattr(self.eventEntry, key, value)

                # Eventdatei einlesen
                with open(str(self.filepfad + '/' + event.attrib['id']), 'r') as csvfile:
                    eventtime = []
                    eventtyp = []
                    for row in csvfile:
                        eventtime.append(row.strip().rsplit(';', 2)[0])
                        eventtyp.append(row.strip().rsplit(';', 2)[1])
                setattr(self.eventEntry, 'event',
                        np.asarray(eventtime, dtype=np.int))
                setattr(self.eventEntry, 'eventtyp',
                        np.asarray(eventtyp, dtype='U10'))

                # EventEntry in Namen der jeweiligen CSV Datei umschreiben
                rename_attribute(self, 'eventEntry',
                                 event.attrib['id'].split('.')[0])

    def getentry(self, name):
        '''Funktion zur Herausgabe von Informationen eines bestimmten Entrytyps'''
        # SignalEntry
        subclassenobject = getattr(self, name)
        if isinstance(subclassenobject, signalEntry):
            return subclassenobject
        # ValuesEntry
        elif isinstance(subclassenobject, valuesEntry):
            return subclassenobject
        # EventEntry
        elif isinstance(subclassenobject, eventEntry):
            return subclassenobject


'''SUBCLASSES'''


class rootAttributes(movisens2python):
    '''RootAttributes'''

    def __init__(self):
        super().__init__
        self.duration = None
        self.measurementId = None
        self.timestampStart = None


class customAttributes(movisens2python):
    '''Custom Attributes'''

    def __init__(self):
        super().__init__
        self.age = None
        self.gender = None
        self.height = None
        self.personId = None
        self.sectorCount = None
        self.sensorLocation = None
        self.sensorSerialNumber = None
        self.sensorType = None
        self.sensorVersion = None
        self.weight = None


class signalEntry(movisens2python):
    '''SignalEntry'''

    def __init__(self):
        super().__init__
        self.adcResolution = None
        self.baseline = None
        self.channels = []
        self.comment = None
        self.contentClass = None
        self.dataType = None
        self.id = None
        self.lsbValue = None
        self.sampleRate = None
        self.unit = None
        self.signal = None


class valuesEntry(movisens2python):
    '''ValuesEntry'''

    def __init__(self):
        super().__init__
        self.adcResolution = None
        self.baseline = None
        self.channels = []
        self.comment = None
        self.contentClass = None
        self.dataType = None
        self.id = None
        self.lsbValue = None
        self.sampleRate = None
        self.unit = None
        self.values = None


class eventEntry(movisens2python):
    '''EventEntry'''

    def __init__(self):
        super().__init__
        self.commentLength = None
        self.id = None
        self.sampleRate = []
        self.typeLength = None
        self.event = None
        self.eventtyp = None


'''HELPERS'''


def rename_attribute(obj, old_name, new_name):
    obj.__dict__[new_name] = obj.__dict__.pop(old_name)


'''MAIN'''


def m2pconverter(*signaltypes, **keywords):
    '''
    Funktion zur Erstellung eines Movisens Objekts

    Non-keyword Arguments
    Welche Signaleinträge gelesen werden sollen:
    *signaltypes = *args [string] [Default = None, entspricht All Data]
    'acc'
    'angularRate'
    'artifact'
    'bpmBxb_live'
    'charging'
    'ecg'
    'hr_live'
    'hrvIsValid_live'
    'hrvRmssd_live'
    'movementAcceleration_live'
    'marker'
    'nn_live'
    'press'
    'stateofCharge'
    'stepCount_live'
    'temp'
    'tempMean_live'

    Keyword Arguments
    **keyowrds = **kwargs [key = string] [Default = None]
    Wenn der XML Baum ausgegeben werden soll:
    showtree = True [Default = False]

    Falls extra weitere Markerdateien/Artifactsdateien genutzt werden:
    extrafile = 'dateiohnepostfix'

    Startzeipunkt der Analyse [ms]: Optional Parameter [string] [Default = None]
    start = 'YYYY-MM-DDThh:mm:ss.xxx'

    Dauer der Anlayse in [s]: Optional Parameter [float] [Default = None]
    dauer = 120123.0
    '''

    # Non-keyword Arguments listen
    datalist = []
    for arg in signaltypes:
        datalist.append(arg)

    # Keyword Arguments listen

    # Times
    if 'start' not in keywords:
        start = None
    else:
        start = keywords['start']
    if 'dauer' not in keywords:
        dauer = None
    else:
        dauer = keywords['dauer']
    if 'extrafile' in keywords:
        extrafile = keywords['extrafile']
    else:
        extrafile = None

    # Objekt erstellen
    movisensobject = movisens2python()
    # Daten auswhählen
    movisensobject.movisens_choosedata()
    # Funktion zum Hinzufügen von Extrafiles
    if extrafile != None:
        movisensobject.movisens_add_content(extrafile)
    # Funktion zur Bearbeitung der Eingabeparameter
    movisensobject.movisens_customsettings(datalist, start, dauer)
    # XML einlesen
    movisensobject.movisens_get_XML()

    # XML Baum zeigen [Default None]
    if 'showtree' in keywords:
        if keywords['showtree'] == True:
            print('Movisensobjekt')
            for key in movisensobject.__dict__:
                if isinstance(movisensobject.__dict__[key], (rootAttributes, customAttributes, signalEntry, valuesEntry, eventEntry)):
                    print(f'-->{key}')
                    for ykey in movisensobject.__dict__[key].__dict__:
                        # if movisensobject.__dict__[key].__dict__[ykey] != None:
                        print(f'-----{ykey}')
                else:
                    # if movisensobject.__dict__[key]:
                    print(f'-----{key}')

    return movisensobject


if __name__ == '__main__':
    
    '''
    EXAMPLE

    (1) Erstellung eines Movisens Objekts mit allen Signalen, die verfügbar sind

    (2) ECG Signalwerte mit LsbValue und Baseline umrechnen (Vorsicht nicht überall vorhanden)

    (3) ECG Signal, Abtastfrequenz, R-Zacken, Marker und Artefakte auslesen

    (4) Signal plotten

    '''

    # Objekt erstellen, mit Signaltyp ECG
    movisensobject = m2pconverter(showtree=True)

    # SignalEntry ECG und ValuesEntry NN_Live auslesen
    ecg = movisensobject.getentry('ecg')
    ecg_signal = (ecg.signal - int(ecg.baseline)) * float(ecg.lsbValue)
    ecg_fs = ecg.sampleRate
    rpeaks = movisensobject.getentry('nn_live').values

    # Plot erstellen

    # ECG
    plt.plot(ecg_signal, label='ECG Signal')

    # R-Zacken
    plt.plot(rpeaks[:, 0], ecg_signal[rpeaks[:, 0]],
             'r+', label='R-Zacken')
   
    # Settings
    plt.title('ECG Single Lead movisens ecgMove4 Chest', fontweight="bold")
    plt.xlabel(f'Samples @ {ecg_fs} ')
    plt.ylabel(f'Amplitude in {ecg.unit} ')

    # Bereiche 10-20 Sekunden
    plt.xlim(10 / (1 / int(ecg_fs)), 20 / (1 / int(ecg_fs)))
    plt.legend()
    plt.show()

