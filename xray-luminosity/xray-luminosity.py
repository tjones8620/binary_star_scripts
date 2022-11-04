import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
plt.style.use('nature')

headers = ['time [s]','L(E>0.1keV)', 'L(E>0.2keV)', 'L(E>0.3keV)', 'L(E>0.5keV)', 'L(E>1keV)', 'L(E>2keV)', 'L(E>5keV)', 'L(E>10keV)']

df = pd.read_csv('/home/visitor_ap4/code/project/scripts/xray-luminosity/xray-lum-wr140-mhd-n256.txt', sep=' ', names=headers, index_col=False)
df['L(2keV < E < 10keV)'] = df['L(E>2keV)'] - df['L(E>10keV)']
df['time [s]'] = df['time [s]'] / (24*60*60)
df.rename(columns={'time [s]': 'time [d]'}, inplace=True)
df.head()

x = df['time [d]']-307.16
# y = df['L(2keV < E < 10keV)']/max(df['L(2keV < E < 10keV)'])

fig, ax = plt.subplots(figsize=(5, 5))
ax.plot(x, df['L(E>0.1keV)'], color='blue', linestyle='solid', label='L(E>0.1keV)')
ax.plot(x, df['L(E>0.2keV)'], color='purple', linestyle='solid', label='L(E>0.2keV)')
ax.plot(x, df['L(E>0.5keV)'], color='red', linestyle='solid', label='L(E>0.5keV)')
ax.plot(x, df['L(E>1keV)'], color='green', linestyle='solid', label='L(E>1keV)')
ax.plot(x, df['L(E>2keV)'], color='orange', linestyle='solid', label='L(E>2keV)')
ax.plot(x, df['L(E>5keV)'], color='cyan', linestyle='solid', label='L(E>5keV)')
ax.plot(x, df['L(E>10keV)'], color='yellow', linestyle='solid', label='L(E>10keV)')

ax.set_xlabel('Days before periastron')
ax.set_ylabel('Unabsorbed Flux [erg/s/cm^2]')
ax.set_title('X-ray Luminosity')
ax.legend()
plt.savefig('/home/visitor_ap4/code/project/scripts/Images/xray-luminosity/xray-lum-wr140-mhd-n256.png', dpi=300)



# def plotting_func(x, y, **kwargs):

#     fig, ax = plt.subplots(figsize = kwargs.get('figsize', (5, 5)))
#     ax.plot(x, y, label=kwargs.get('label', ''))
#     ax.set_xlabel(kwargs.get('xlabel', ''))
#     ax.set_ylabel(kwargs.get('ylabel', ''))
#     ax.set_title(kwargs.get('title', ''))
#     ax.set_xlim(kwargs.get('xlim', None))
#     ax.set_ylim(kwargs.get('ylim', None))


#     if kwargs.get('save', True):
#         plt.savefig(os.path.join(kwargs.get('img_dir', os.path.abspath(os.path.dirname(__file__))), kwargs.get('img_name', f"image{np.random.randint(1000, size=1)}")), dpi=kwargs.get('dpi', 300))

#     if kwargs.get('show', True):
#         plt.show()
    
#     return fig, ax

# plotargs = {'xlabel': 'Time [d]', 'ylabel': 'L(2keV < E < 10keV)', 
#             'save': True, 'title': 'X-ray Luminosity', 'img_name': 'xray-lum-wr140-mhd-n256.png', 
#             'img_dir': '/home/visitor_ap4/code/project/scripts/Images/xray-luminosity',
#             'xlim': (0, 100), 'ylim': (0, 1.1)
#             }

# energy_headers = ['L(E>0.1keV)', 'L(E>0.2keV)', 'L(E>0.3keV)', 'L(E>0.5keV)', 'L(E>1keV)', 'L(E>2keV)', 'L(E>5keV)', 'L(E>10keV)']


