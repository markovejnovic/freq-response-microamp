#!/usr/bin/python
from oscilloscope_waveform import OscilloscopeWaveform
import matplotlib.pyplot as plt
import numpy as np
import sys
from glob import glob
import os
import csv
import tty
import termios
import re
import time

def input_char(message):
    return input(message)

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [atoi(c) for c in re.split('(\d+)', text)]

if __name__ == '__main__':
    recursive = False
    verbose = False

    if '-d' in sys.argv or '--debug' in sys.argv or '-v' in sys.argv or \
            '--verbose' in sys.argv:
        verbose = True
    if '-r' in sys.argv or '--recursive' in sys.argv:
        recursive = True
    directory = sys.argv[-1]

    if recursive == False:
        print('Non recursive behavior is not supported. Yeah...')
    elif recursive == True:
        with open('data.csv', 'w', newline='') as data_csv:
            csvwriter = csv.writer(data_csv, delimiter=',')
            csvwriter.writerow(['Scan Name', 'Source', 'Time', 'Voltage', 
                'Fit Amplitude', 'Fit Frequency', 'Fit Phase Shift',
                'Fit Offset'])
            np.set_printoptions(threshold=np.Inf)
            for filename in sorted(glob(os.path.join(directory, '**', '*.CSV'), 
                    recursive=True), key=natural_keys):
                print('[INFO] Analyzing ' + filename + '...')

                ow = OscilloscopeWaveform(filename)
                guesses = ow.guess_fit_values()
                if guesses[4] == 0:
                    to_write = False
                else:
                    to_write = True
                ow.best_fit(auto_guess=False, 
                        guess_mean=guesses[0],
                        guess_std=guesses[1],
                        guess_phase=guesses[2],
                        guess_freq=guesses[3],
                        guess_amp=guesses[4], 
                        force=True)

                if to_write:
                    csvwriter.writerow([filename, ow.source, \
                            re.sub(' +', ' ', 
                                np.array_repr(ow.x).replace('\n', '')), \
                            re.sub(' +', ' ', 
                                np.array_repr(ow.y).replace('\n', '')), \
                            ow.best_fit()[0], ow.best_fit()[1], \
                            ow.best_fit()[2], ow.best_fit()[3]])
                print('[INFO] FINAL FIT')
                print('Frequency: ' + str(ow.best_fit()[1]))
                print('Amplitude: ' + str(ow.best_fit()[0]) + '\n')
                data_csv.flush()
                os.fsync(data_csv.fileno())
                print('[INFO] Finished ' + filename + '.\n')
