import wfdb
import numpy as np
import os
from tkinter import filedialog, Tk


class phdb2python():
    '''
    Klasse zum einlesen und bereitstellen von PHDB-Datensätzen
    '''
    # Verbesserungspotential:
    # - Generalisieren: Datenformate der WDFB Datenbank
    # - Anstatt ADC gain, ADC res.. Eigenes Dictionary Einführen
    # EXAMPLE SCHREIBEN wie mov2python

    def __init__(self):
        ''' Datenobjekt PHDB Daten '''
        self.filename = None
        self.filepfad = None
        self.name = None
        self.signal = []
        self.marker = []
        self.fs = 1
        self.adc_gain = None
        self.adc_res = None
        self.adc_zero = None
        self.units = None

    def choosedata(self, extension=['*.atr;*.ari'], datatype='PHDB Dataset Files'):
        '''Funktion zum Auswählen der Daten'''

        root = Tk()  # TinkerDialog öffnen
        root.withdraw()  # Tinkerfenster verstecken
        files = filedialog.askopenfilenames(
            defaultextension='.*',
            filetypes=[(datatype, extension), ],
            title="Select Measurement",
            multiple=True,
        )
        self.filepfad = os.path.dirname(files[0])
        if len(files) <= 1:
            self.filename = os.path.split(files[0])[1]
        else:
            self.filename = os.path.split(self.filepfad)[1]
        self.name = os.path.splitext(self.filename)[0]

    def phdb_read_ecg_channel(self):
        '''Funktion zum Einlesen eines ECG Channel von PHDB Daten '''
        self.signal = wfdb.rdrecord(
            str(self.filepfad + '/' + self.name), channels=[0])
        attitudes = ['adc_gain', 'adc_res', 'adc_zero', 'fs', 'units']
        for attr in attitudes:
            if hasattr(self.signal, attr):
                setattr(self, attr, getattr(self.signal, attr))

    def phdb_read_marker(self):
        '''Funktion zum Einlesen von Markern der PHDB Daten '''
        annotationtype = self.filename.split('.')[1]
        self.marker = wfdb.rdann(
            str(self.filepfad) + '/' + self.name, annotationtype)

    def phdb_make_ecg_great_again(self):
        '''Funktion zum Umwandeln eines ECG Signals in einen Objectinternen Numpy Signalarray '''
        signalwert = []
        for line in self.signal.p_signal:
            signalwert.append(line[0])
        self.signal = np.asarray(signalwert)

    def phdb_beatmarker_only(self):
        '''
        Funktion zum Umwandeln der Markerinformationen in einen Objectinternen Markerarray
        - Herausgefilter werden Marker, die keinen Herzschlag markieren
         '''
        phdb_notbeatlist = ['[', '!', ']', 'x', '(', ')', 'p', 't', 'u', '`',
                            "'", '^', '|', '~', '+', 's', 'T', '*', 'D', '=', '"', '@']
        newmarker = []
        for idx, k in enumerate(self.marker.sample):
            if self.marker.symbol[idx] in phdb_notbeatlist:
                continue
            else:
                newmarker.append(k)
        self.marker = np.asarray(newmarker)


def phdb2pconverter(showtree=False):
    phdbobject = phdb2python()
    phdbobject.choosedata()
    # ECG Einlesen
    phdbobject.phdb_read_ecg_channel()
    # Annotationen einlesen
    phdbobject.phdb_read_marker()
    # ECG Signal umwandeln
    phdbobject.phdb_make_ecg_great_again()
    # Markerfilter: R-Zacken
    phdbobject.phdb_beatmarker_only()

    if showtree == True:
        print(phdbobject.__dict__)

    return phdbobject


if __name__ == "__main__":
    phdbdata = phdb2pconverter(showtree=True)
