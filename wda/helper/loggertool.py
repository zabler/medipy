import time
import sys
import os


def initlog(logfile):
    logdatei = str(os.path.abspath(".")) + '/log/' + str(logfile) + '.log '
    if not os.path.exists(str(os.path.abspath(".")) + '/log'):
        os.makedirs(str(os.path.abspath(".")) + '/log')
    if not os.path.exists(logdatei):
        newfile = open(logdatei, 'w')
        newfile.write("%s %s\n" % (time.strftime(
            "%d.%m.%Y %H:%M:%S"), 'Init Logger ' + str(logfile)))
        newfile.close
    return logdatei


def log(msg, logfile):
    '''Save Messages as Logdatei
    Input Message und Log-Datei
    Return Print Message und Speichern in Logdatei '''
    print(msg)
    msg = str(msg)
    with open(logfile, 'a') as logger:
        logger.write("%s %s\n" % (time.strftime("%d.%m.%Y %H:%M:%S"), msg))

if __name__ == "__main__":
    
    print('''
    No Example: Funktion in Script einbinden um Infomartionen in Logdatei zu speichern
    ''')

