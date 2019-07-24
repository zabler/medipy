from wfdb import processing as wfdbprocessor

''' Funktionen zum Auswertung von Human Reviewed vs Autodetected'''
# Verbesserung
# - eigene auswettungsclasseVergleichsfunktion, Plotfunktion und Stastik aufbauen TP, TN , ..


def mitcompare(signal, refdetections, calcdetections, sampleaccuracy=50):
    '''Funktion zum Auswerten der Daten mittels Statistikfunktionen des WFDB Packets '''
    # Vergleich der Algorithmen mittels Compare Klasse des WFDB Tools
    compartitor = wfdbprocessor.compare_annotations(
        refdetections, calcdetections, sampleaccuracy, signal)
    compartitor.print_summary()
    compartitor.plot()


def owncompaore(signal, refdetections, calcdetections, sampleaccuracy):

    pass


def ownstats(TP, TN, FP, FN):
    '''Stats
    Tabelle Referenz: Positive Negative
    Detektion Positive   TP     FP  
    Detektion Negative   FN     TN

    TP True positive: Detektionen bei denen ein Referenzwert vorliegt
    FP False positive: Detektionen bei denen kein Referenzwert vorliegt
    FN False negative: Keine Detektionen bei denen ein Referenzwert vorliegt
    TN True Negative: Keine Detektionen be denen kein Referenzwert vorliegt

    N_dec = Anzahl der Detektionen = TP + FP
    N_ref = Anzahl der Referenzen = TP + FN

    Spezifizität = TP/N_ref = TP/(TP+FN)
    Positive Predicitivity = TP/N_dec = TP/(TP+FP)
    Falschalarmrate = FP/N_dec = FP/(TP+FP)
    Genauigkeit = TP/(TP+FP+FN)


    '''
    pass
