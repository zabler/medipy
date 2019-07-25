from wfdb import processing as wfdbprocessor

def mitcompare(signal, refdetections, calcdetections, sampleaccuracy=50):
    '''Funktion zum Auswerten der Daten mittels Statistikfunktionen des WFDB Packets '''
    # Vergleich der Algorithmen mittels Compare Klasse des WFDB Tools
    compartitor = wfdbprocessor.compare_annotations(
        refdetections, calcdetections, sampleaccuracy, signal)
    compartitor.print_summary()
    compartitor.plot()


if __name__ == "__main__":
    
    print('''
    No Example: Funktion in Script einbinden um die Genauigkeit des Algorithmus zu prüfen
    ''')

