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
            plt.ion()
            csvwriter = csv.writer(data_csv, delimiter=',')
            csvwriter.writerow(['Scan Name', 'Source', 'Time', 'Voltage', 
                'Fit Amplitude', 'Fit Frequency', 'Fit Phase Shift',
                'Fit Offset'])
            freq = 100
            amp = 2.5
            for filename in sorted(glob(os.path.join(directory, '**', '*.CSV'), 
                    recursive=True), key=natural_keys):
                print('[INFO] Analyzing ' + filename + '...')
                ow = OscilloscopeWaveform(filename)
                ow.best_fit(guess_freq=freq, guess_amp=amp, force=True)
                ow.plot(plot_fit=True)
                plt.show()
                c = input_char("[INFO] Press A (accept) or Z (deny) to " +
                        "continue...")
                #c = 'a'
                if c == 'Z' or c == 'z':
                    while(c is 'Z' or c is 'z'):
                        freq = float(input('[INFO] Input frequency guess: '))
                        amp = 2.8
                        plt.close()
                        ow.best_fit(guess_freq=freq, guess_amp=amp, 
                                auto_guess=False, force=True)
                        ow.plot(plot_fit=True)
                        plt.show()
                        c = input_char('[INFO] Press A (accept) or Z (deny)' + 
                                ' continue...')
                else:
                    np.set_printoptions(threshold=np.Inf)
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
                plt.close()
                print('[INFO] Finished ' + filename + '.\n')
