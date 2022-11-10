import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
plt.style.use("ieee")
import pathlib
home = pathlib.Path.home()

############################################################################################################

headers = ['time [s]','L(E>0.1keV)', 'L(E>0.2keV)', 'L(E>0.3keV)', 'L(E>0.5keV)', 'L(E>1keV)', 'L(E>2keV)', 'L(E>5keV)', 'L(E>10keV)']

path = os.path.join(home, "project/scripts/xray-luminosity/xray-lum-wr140-mhd-n256.txt")
df = pd.read_csv(path, sep=' ', names=headers, index_col=False)

df['L(2keV < E < 10keV)'] = df['L(E>2keV)'] - df['L(E>10keV)']
df['time [s]'] = df['time [s]'] / (24*60*60)
df.rename(columns={'time [s]': 'time [d]'}, inplace=True)

############################################################################################################

# Plotting xray luminosity vs time for different energy bands
fig, ax1 = plt.subplots(figsize=(6,4))

plt.title("Unabsorbed Flux vs Time")
ax1.set_xlabel('Time (d)')
ax1.set_ylabel('Unabsorbed Flux ($erg$ $cm^{-2}$$s^{-1}$)')

ax1.grid(True)

ax2 = ax1.twinx()  # create a second axes that shares the same x-axis

ax2.plot(df['time [d]']-307.16, df['L(E>0.1keV)'], color='blue', linestyle='solid', label='$L(E>0.1keV)$')
ax2.plot(df['time [d]']-307.16, df['L(E>0.2keV)'], color='purple', linestyle='solid', label='$L(E>0.2keV)$')
ax2.plot(df['time [d]']-307.16, df['L(E>0.5keV)'], color='red', linestyle='solid', label='$L(E>0.5keV)$')
ax2.plot(df['time [d]']-307.16, df['L(E>1keV)'], color='yellow', linestyle='solid', label='$L(E>1keV)$')
ax2.plot(df['time [d]']-307.16, df['L(E>2keV)'], color='cyan', linestyle='solid', label='$L(E>2keV)$')
ax2.plot(df['time [d]']-307.16, df['L(E>5keV)'], color='lightpink', linestyle='solid', label='$L(E>5keV)$')
ax2.plot(df['time [d]']-307.16, df['L(E>10keV)'], color='peru', linestyle='solid', label='$L(E>10keV)$')

ax2.set_ylabel('X-ray Luminosity ($erg$ $s^{-1}$)')
ax2.legend(loc='upper left', ncol=2, framealpha=1, frameon=True, fontsize=5)
plt.savefig('xray-lum-wr140-mhd-n256.png', dpi=300, bbox_inches='tight')

############################################################################################################

# Plotting xray luminosity vs time for 2-10keV band

fig, ax1 = plt.subplots(figsize=(6,4))

plt.title("Unabsorbed Flux vs Time")
ax1.set_xlabel('Time (d)')
ax1.set_ylabel('Unabsorbed Flux ($erg$ $cm^{-2}$$s^{-1}$)')

ax1.grid(True)

ax2 = ax1.twinx()  # create a second axes that shares the same x-axis

daysec = 24 * 60 * 60 

x = df['time [d]']-307.16
y = df['L(2keV < E < 10keV)'] #/max(df['L(2keV < E < 10keV)'])

ax2.plot(x, y, color='blue', linestyle='solid', label='$L(2keV - 10keV)$')

ax2.set_ylabel('X-ray Luminosity ($erg$ $s^{-1}$)')
ax2.legend(loc='upper left', ncol=2, framealpha=1, frameon=True, fontsize=5)
plt.savefig('xray-lum-wr140-mhd-n256-2_10keV.png', dpi=300, bbox_inches='tight')