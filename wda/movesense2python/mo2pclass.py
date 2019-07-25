import numpy as np
from tkinter import filedialog, Tk
import os
import csv


class movesense2python():
    def __init__(self):
        '''Init der BaseClass Objektvariablen'''
        # Choosedata
        self.filename = None
        self.filepfad = None
        self.name = None

    def choosedata(self, extension='.csv', datatype='Movesense File'):
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

    def ecgreader(self):
        with open(str(self.filepfad + '/' + self.filename), 'r') as csvfile:
            read = list(csv.reader(csvfile, delimiter=','))
            ecgtime = []
            ecg = []
            for k in read:
                ecgtime.append(k[0])
                ecg.append(k[1])
            self.ecgtime = np.array(ecgtime[1:], dtype=np.int)
            self.ecg = np.array(ecg[1:], dtype=np.int)

    def rrreader(self):
        with open(str(self.filepfad + '/' + self.filename), 'r') as csvfile:
            read = list(csv.reader(csvfile, delimiter=','))
            # self.rr = np.array(read[1:])  # , dtype=np.float16)
            rtime = []
            hr = []
            rr = []
            for k in read:
                rtime.append(k[0])
                hr.append(k[1])
                rr.append(k[2])
            self.rtime = np.array(rtime[1:], dtype=np.int)
            self.hr = np.array(hr[1:], dtype=np.float16)
            self.rr = np.array(rr[1:], dtype=np.int)

    '''Further Data'''
    # def imu9reader(self):
    #     pass


def mo2pconverter(*signaltypes):
    # Objekt erstellen
    movesenseobject = movesense2python()
    # Daten auswhählen
    movesenseobject.choosedata()
    # Welche Daten werden eingelesen
    if 'heartrate' in signaltypes:
        movesenseobject.rrreader()
    elif 'ecg' in signaltypes:
        movesenseobject.ecgreader()
    return movesenseobject


if __name__ == '__main__':
     '''
    EXAMPLE

    (1) Movesense Object erstellen

    (2) Movesense HR Datei auswählen

    (3) RR Reader Movesense liest HR Daten ein

    (4) Print der RR-Intervalle

    '''
    #Movesense Object erstellen
    movesenseobject = movesense2python()

    #Movesense HR Datei auswählen
    movesenseobject.choosedata()

    #RR Reader Movesense liest HR Daten ein
    movesenseobject.rrreader()
    
    #Print der RR-Intervalle
    print(movesenseobject.rr)
