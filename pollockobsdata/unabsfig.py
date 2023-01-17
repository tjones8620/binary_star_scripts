from matplotlib.pyplot import *
import pandas as pd
import astropy.units as u
from astropy.time import Time
import astropy.constants as c
import numpy as np
import os

Period = 2895.00*u.day
To = 60636.23
epoch = Time(To-5*Period.value, format='mjd')
eccentricity = 0.8993
omega=227.44
m1 = 10.31 
m2 = 29.27
mtot=(m1+m2)*u.Msun 
sma = (((c.G*mtot.to(u.kg)/(4*np.pi**2)*(Period)**2).decompose())**(1./3)).to(u.AU)
tacalc = np.arange(201)*2*np.pi/200
sepcalc = sma*(1-eccentricity**2)/(1+eccentricity*np.cos(tacalc))
D = 1518*u.pc

cwd = os.path.dirname(os.path.abspath(__file__))


rxtexlnamecor = 'Fix_excel_spreadsheet_20200916_2.0_10.0_rxte.xls'
swiftxlnamecor = 'Fix_excel_spreadsheet_20200916_2.0_10.0_swift.xls'
nicerxlnamecor = 'Fix_excel_spreadsheet_20200916_2.0_10.0_nicer.xls'

rxteunabs = pd.read_excel(os.path.join(cwd, rxtexlnamecor), index_col=0)
swiftunabs = pd.read_excel(os.path.join(cwd,swiftxlnamecor), index_col=0)
nicerunabs = pd.read_excel(os.path.join(cwd, nicerxlnamecor), index_col=0)

fig, ax = subplots(figsize=[15,8])

x1= -3500
x2 = 4500
xlimit = [x1, x2]
ylimit = [0,3]

ax.set_xlabel('Time in Days relative to the Periastron passage of {0}'.format('2009-01-15'),fontsize=16)
ax.set_ylabel('Unabsorbed L$_{x}$ (10$^{34}$ ergs s$^{-1}$)',fontsize=16)
ax.set_ylim(ylimit)
ax.set_xlim(xlimit)

scl=1e11
jdoff = (epoch+3*Period).jd

i=0

minexpo = 200
lumscl=4*np.pi*D.to(u.cm)**2/1e34
rbkg =0.5
wtbkg = .8 # for WT mode data
pcbkg = 0.0
nbkg = 0.7


ax.plot(rxteunabs.JDMID-jdoff, (rxteunabs.UnabsFlux/1.35 - rbkg/scl)*lumscl,'ko', label='RXTE')
t_rxte = rxteunabs.JDMID-jdoff
l_rxte = (rxteunabs.UnabsFlux/1.35 - rbkg/scl)*lumscl

pc = (swiftunabs['mode'] == 'pc')
ax.plot(swiftunabs[pc].JDMID-jdoff, swiftunabs[pc].UnabsFlux*lumscl,'bo',label='Swift PC')
t_swiftpc = swiftunabs[pc].JDMID-jdoff
l_swiftpc = swiftunabs[pc].UnabsFlux*lumscl

pc = (swiftunabs['mode'] == 'wt')
ax.plot(swiftunabs[pc].JDMID-jdoff, (swiftunabs[pc].UnabsFlux-wtbkg/scl)*lumscl,'bo',label='Swift WT', markerfacecolor='none')
t_swiftwt = swiftunabs[pc].JDMID-jdoff
l_swiftwt = (swiftunabs[pc].UnabsFlux-wtbkg/scl)*lumscl

ax.plot(nicerunabs.JDMID-jdoff, (nicerunabs.UnabsFlux - nbkg/scl)*lumscl,'go', label='NICER')
t_nicer = nicerunabs.JDMID-jdoff
l_nicer = (nicerunabs.UnabsFlux - nbkg/scl)*lumscl

# Create a dataframe with all the data
df = pd.DataFrame({'t_rxte':t_rxte, 'l_rxte':l_rxte, 't_swiftpc':t_swiftpc, 'l_swiftpc':l_swiftpc, 't_swiftwt':t_swiftwt, 'l_swiftwt':l_swiftwt, 't_nicer':t_nicer, 'l_nicer':l_nicer})
print(df.head(50))
df.to_csv('wr140unabslum.csv')

# sepnorm1=7

# ax.plot(np.asarray(pcalc)*Period.value, 1./sepcalc*sepnorm1, color='grey', linestyle='--', label='1/D')
# ax.plot(np.asarray(pcalc)*Period.value-Period.value, 1./sepcalc*sepnorm1,  color='grey', linestyle='--')
# ax.plot(np.asarray(pcalc)*Period.value+Period.value, 1./sepcalc*sepnorm1,  color='grey', linestyle='--')
# ax.plot(np.asarray(pcalc)*Period.value-2*Period.value, 1./sepcalc*sepnorm1,  color='grey', linestyle='--')

legend(loc='upper left', ncol = 2, fontsize=12, numpoints=1, frameon=True, facecolor='white', framealpha=1)

# Add phase axis on Top
# topax = ax.twiny()
# p1 = x1/Period.value
# p2 = x2/Period.value
# topax.set_xlim([p1, p2])
# topax.set_xlabel('Phase', fontsize=20)
# topax.tick_params(axis='x', labelsize=16 )
# plot()

plotdir = "./"

save=False
if save:
    figname=os.path.join(plotdir,'unabsorbed_luminosity.png')
    print(f'Writing {figname}')
    savefig(figname)