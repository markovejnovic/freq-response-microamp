#!/usr/bin/python
import csv
import sys
import matplotlib.pyplot as plt
from math import pi
from numpy import log10

freqs = []
amp = []
const = []

f = sys.argv[-1]

with open(f, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    skip_line = True
    for row in reader:
        if skip_line:
            skip_line = False
            continue
        freqs.append(float(row[1]))
        amp.append(float(row[0]))

freqs, amp = (list(t) for t in zip(*sorted(zip(freqs, amp))))

freqs_final = []
for f in freqs:
    freqs_final.append(f / 2 / pi)

amp_final = []
for a in amp:
    amp_final.append(20 * log10(abs(a)/0.15))

for f in freqs:
    const.append(0.15)

plt.grid()
plt.xlabel('Frequency [Hz]')
plt.ylabel('Gain [dB]')
plt.plot(freqs_final, amp_final, 'r')
plt.xscale('log')
axes = plt.gca()
plt.show()
