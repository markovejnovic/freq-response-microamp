import csv
from enum import Enum
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import leastsq

class OscilloscopeWaveform:
    def __init__(self, data, memory_length, trigger_level, source, probe, 
            vertical_units, vertical_scale, horizontal_units, 
            horizontal_scale, horizontal_position, horizontal_mode, 
            sampling_period, firmware, time, mode, osc_vsampling_max=100, 
            osc_vsections=10):
        self.data = data
        self.memory_length = memory_length
        self.trigger_level = trigger_level
        self.source = source
        self.probe = probe
        self.vertical_units = vertical_units
        self.vertical_scale = vertical_scale
        self.horizontal_units = horizontal_units
        self.horizontal_scale = horizontal_scale
        self.horizontal_position = horizontal_position
        self.horizontal_mode = horizontal_mode
        self.sampling_period = sampling_period
        self.firmware = firmware
        self.time = time
        self.mode = mode
        self.osc_vsampling_max = osc_vsampling_max
        self.vsections = vsections

        self._best_fit = None

    def __init__(self, f, osc_vsampling_max=100, osc_vsections=10):
        """Overloaded constructor which allows for reading from a .csv file so 
        this class automatically figures its stuff out. Look, I'm not writing 
        this for a living.
        """
        self.osc_vsampling_max = osc_vsampling_max
        self.osc_vsections = osc_vsections
        with open(f, newline='') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            analyzing_data = False # Assume the file starts with the header
            self.data = []
            self._best_fit = None
            for row in csvreader:
                if analyzing_data == True:
                    self.data = np.append(self.data, float(row[0]))
                else:
                    if row[0] == 'Memory Length':
                        self.memory_length = float(row[1])
                    elif row[0] == 'Trigger Level':
                        self.trigger_level = float(row[1])
                    elif row[0] == 'Source':
                        self.source = row[1]
                    elif row[0] == 'Probe':
                        self.probe = float(row[1][:-1])
                    elif row[0] == 'Vertical Units':
                        self.vertical_units = row[1]
                    elif row[0] == 'Vertical Scale':
                        self.vertical_scale = float(row[1])
                    elif row[0] == 'Vertical Position':
                        self.vertical_position = float(row[1])
                    elif row[0] == 'Horizontal Units':
                        self.horizontal_units = row[1]
                    elif row[0] == 'Horizontal Scale':
                        self.horizontal_scale = float(row[1])
                    elif row[0] == 'Horizontal Position':
                        self.horizontal_position = float(row[1])
                    elif row[0] == 'Horizontal Mode':
                        self.horizontal_mode = row[1]
                    elif row[0] == 'Sampling Period':
                        self.sampling_period = float(row[1])
                    elif row[0] == 'Firmware':
                        self.firmware = float(row[1][1:])
                    elif row[0] == 'Time':
                        self.time = row[1]
                    elif row[0] == 'Mode':
                        self.mode = row[1]
                    elif row[0] == 'Waveform Data':
                        analyzing_data = True
    
    @property
    def data(self):
        """The waveform data"""
        return self._data

    @data.setter
    def data(self, data):
        if isinstance(data, np.ndarray):
            self._data = data
        elif isinstance(data, list):
            self._data = np.array(data)
        else:
            raise Exception('Incorrect data type. Please provide a np.array.')

    @property
    def memory_length(self):
        """The memory length"""
        return self._memory_length

    @memory_length.setter
    def memory_length(self, memory_length):
        self._memory_length = memory_length

    @property
    def trigger_level(self):
        """Description of trigger_level"""
        return self._trigger_level

    @trigger_level.setter
    def trigger_level(self, trigger_level):
        self._trigger_level = trigger_level

    @property
    def source(self):
        """Description of source"""
        return self._source

    @source.setter
    def source(self, source):
        self._source = source
    
    @property
    def probe(self):
        """Description of probe"""
        return self._probe

    @probe.setter
    def probe(self, probe):
        self._probe = probe

    @property
    def vertical_units(self):
        """Description of vertical_units"""
        return self._vertical_units

    @vertical_units.setter
    def vertical_units(self, vertical_units):
        self._vertical_units = vertical_units

    @property
    def vertical_scale(self):
        """Description of vertical_scale"""
        return self._vertical_scale

    @vertical_scale.setter
    def vertical_scale(self, vertical_scale):
        self._vertical_scale = vertical_scale

    @property
    def vertical_position(self):
        """Description of vertical_position"""
        return self._vertical_position

    @vertical_position.setter
    def vertical_position(self, vertical_position):
        self._vertical_position = vertical_position

    @property
    def horizontal_units(self):
        """Description of horizontal_units"""
        return self._horizontal_units

    @horizontal_units.setter
    def horizontal_units(self, horizontal_units):
        self._horizontal_units = horizontal_units

    @property
    def horizontal_scale(self):
        """Description of horizontal_scale"""
        return self._horizontal_scale

    @horizontal_scale.setter
    def horizontal_scale(self, horizontal_scale):
        self._horizontal_scale = horizontal_scale

    @property
    def horizontal_position(self):
        """Description of horizontal_position"""
        return self._horizontal_position

    @horizontal_position.setter
    def horizontal_position(self, horizontal_position):
        self._horizontal_position = horizontal_position

    @property
    def horizontal_mode(self):
        """Description of horizontal_mode"""
        return self._horizontal_mode

    @horizontal_mode.setter
    def horizontal_mode(self, horizontal_mode):
        self._horizontal_mode = horizontal_mode

    @property
    def sampling_period(self):
        """Description of sampling_period"""
        return self._sampling_period

    @sampling_period.setter
    def sampling_period(self, sampling_period):
        self._sampling_period = sampling_period

    @property
    def firmware(self):
        """Description of firmware"""
        return self._firmware

    @firmware.setter
    def firmware(self, firmware):
        self._firmware = firmware

    @property
    def time(self):
        """Description of time"""
        return self._time

    @time.setter
    def time(self, time):
        self._time = time

    @property
    def mode(self):
        """Description of mode"""
        return self._mode

    @mode.setter
    def mode(self, mode):
        self._mode = mode

    @property
    def osc_vsampling_max(self):
        """The maximum sampling point of the oscilloscope"""
        return self._osc_vsampling_max

    @osc_vsampling_max.setter
    def osc_vsampling_max(self, osc_vsampling_max):
        self._osc_vsampling_max = osc_vsampling_max

    @property
    def osc_vsections(self):
        """The number of vertical sections the oscilloscope has"""
        return self._osc_vsections

    @osc_vsections.setter
    def osc_vsections(self, osc_vsections):
        self._osc_vsections = osc_vsections

    @property
    def figure(self):
        """A pyplot figure"""
        return self._figure

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    def calculate_x_y(self):
        self._x = np.linspace(0, self.data.size * self.sampling_period, 
                self.data.size)
        self._y = (self.data * self.vertical_scale + self.vertical_position) / \
                self.osc_vsampling_max * self.osc_vsections / 2

    def plot(self, plot_fit=False, ylim=[-5, 5]):
        """Plots a pyplot.

        Does not call matplotlib.pyplot.show(), you are required to do this
        """
        self.calculate_x_y()

        self._figure, ax = plt.subplots()
        ax.plot(self.x, self.y)
        ax.set_xlabel('Time [' + self.horizontal_units + ']')
        ax.set_ylabel('Voltage [' + self.vertical_units + ']')
        ax.set_ylim(ylim)

        if plot_fit:
            fit = self.best_fit()
            ax.plot(self.x, fit[0] * np.sin(fit[1] * self.x + fit[2]) + \
                    fit[3])

    def _peak_det(self, v=None, delta=0.5, x = None):
        maxtab = []
        mintab = []

        self.calculate_x_y()

        if v is None:
            v = self._y
           
        if x is None:
            x = self._x
        
        v = np.asarray(v)
        
        if len(v) != len(x):
            raise Exception('Input vectors v and x must have same length')
        
        if not np.isscalar(delta):
            raise Exception('Input argument delta must be a scalar')
        
        if delta <= 0:
            raise Exception('Input argument delta must be positive')
        
        mn, mx = np.Inf, -np.Inf
        mnpos, mxpos = np.NaN, np.NaN
        
        lookformax = True
        
        for i in np.arange(len(v)):
            this = v[i]
            if this > mx:
                mx = this
                mxpos = x[i]
            if this < mn:
                mn = this
                mnpos = x[i]
            
            if lookformax:
                if this < mx-delta:
                    maxtab.append((mxpos, mx))
                    mn = this
                    mnpos = x[i]
                    lookformax = False
            else:
                if this > mn+delta:
                    mintab.append((mnpos, mn))
                    mx = this
                    mxpos = x[i]
                    lookformax = True

        return np.array(maxtab), np.array(mintab)

    def guess_fit_values(self):
        """Tries to roughly guess fit values for the sine wave"""
        self.calculate_x_y()

        guess_mean = np.mean(self.y)
        guess_std = 3 * np.std(self.y) / (2**0.5) / (2**0.5)
        guess_phase = 0
        max_peaks, min_peaks = self._peak_det()
        if not len(max_peaks):
            guess_freq = 0
        else:
            second = max_peaks[1][0]
            second_to_last = max_peaks[-1][0]
            peaks_n = max_peaks.T[0].size - 2
            guess_freq = (1 / ((second_to_last - second)/peaks_n)* 2 * np.pi) 
        print('Guessed frequency: ' + str(guess_freq))
            
        guess_amp = np.amax(self.y)
        return [guess_mean, guess_std, guess_phase, guess_freq, guess_amp]

    def best_fit(self, guess_mean=None, guess_std=None, guess_phase=0,
            guess_freq=1200, guess_amp=3, force=False, auto_guess=True):
        """Calculates the best_fit value if it is not calculated, otherwise 
        returns the already calculated one

        Returns:
            list - the best fit coefficients in 0 * sin(1 * x + 2) + 3 form
        """
        self.calculate_x_y()
        if guess_mean == None:
            guess_mean = np.mean(self.y)
        if guess_std == None:
            guess_std = 3 * np.std(self.y) / (2**0.5) / (2**0.5)

        if auto_guess:
            guesses = self.guess_fit_values()
            guess_mean = guesses[0]
            guess_std = guesses[1]
            guess_phase = guesses[2]
            guess_freq = guesses[3]
            guess_amp = guesses[4]

        if self._best_fit is None or force == True:
            optimize_func = lambda x: x[0] * np.sin(x[1] * self.x + x[2]) + \
                    x[3] - self.y
            est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, 
                    [guess_amp, guess_freq, guess_phase, guess_mean])[0]
            self._best_fit = [est_amp, est_freq, est_phase, est_mean]
        return self._best_fit
