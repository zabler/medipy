import os
import wfdb
from tkinter import filedialog, Tk


class phdbloader():

    def __init__(self):
        '''Dataobjekt Parameter'''
        self.filepfad = None
        self.db = None

    def choosedatabase(self):
        # Datenbanken zeigen
        print(wfdb.io.get_dbs())
        # Welche Datenbank?
        print('_________________________________________________________')
        self.db = input('Welche Datenbank soll verwendet werden? xxxxx: ')

    def choosesavepath(self):
        '''Funktion zum Auswählen eines einzelnen Dateipfades'''
        # Wo soll gespeichert werden
        root = Tk()
        root.withdraw()
        self.filepfad = filedialog.askdirectory(
            title="Speicherort wählen",)
        self.name = os.path.dirname(self.filepfad)

    def phdb_download(self):
        '''Funktion zum Download von PHDB Daten '''
        wfdb.dl_database(self.db + '/', self.filepfad + '/')


def databaseloader():

    download = phdbloader()
    download.choosedatabase()
    download.choosesavepath()
    download.phdb_download()


if __name__ == '__main__':
    databaseloader()
