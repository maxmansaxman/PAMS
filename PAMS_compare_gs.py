'''Compares methane clumped values determined using adduct lines and an in situ goal seek method'''
import matplotlib as mpl
# mpl.use('pdf')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

a = pd.read_csv('autoFlatListExport_august2016_adduct.csv', dialect = 'excel')
a_gs = pd.read_csv('autoFlatListExport_august2016_gs.csv', dialect = 'excel')

#dD comparison
fig, ax = plt.subplots()
ax.errorbar(a['d2H (vsmow)'],a_gs['d2H (vsmow)'] - a['d2H (vsmow)'], xerr= a['d2H_sterr'].values*2, yerr = np.sqrt(a_gs['d2H_sterr']**2+a['d2H_sterr']**2).values*2, fmt = 'bs' )
ax.plot(np.linspace(a['d2H (vsmow)'].min(), a['d2H (vsmow)'].max()), np.zeros(50), 'g-')
ax.set_xlabel(ur'$\delta\mathrm{^{2}H_{vsmow, adductLine}}$ (\u2030)')
ax.set_ylabel(ur'$\Delta\mathrm{^{2}H_{goalSeek - adductLine}}$ (\u2030)')
plt.savefig('d2H_gs_adduct_comparison.pdf')

fig, ax = plt.subplots()
ax.errorbar(a['d13C (vpdb)'],a_gs['d13C (vpdb)'] - a['d13C (vpdb)'], xerr= a['d13C_sterr'].values*2, yerr = np.sqrt(a_gs['d13C_sterr']**2+a['d13C_sterr']**2).values*2, fmt = 'bs' )
ax.plot(np.linspace(a['d13C (vpdb)'].min(), a['d13C (vpdb)'].max()), np.zeros(50), 'g-')
ax.set_xlabel(ur'$\delta\mathrm{^{13}C_{vpdb, adductLine}}$ (\u2030)')
ax.set_ylabel(ur'$\Delta\mathrm{^{13}C_{goalSeek - adductLine}}$ (\u2030)')
plt.savefig('d13C_gs_adduct_comparison.pdf')

fig, ax = plt.subplots()
ax.errorbar(a['D18_stochastic'],a_gs['D18_stochastic'] - a['D18_stochastic'], xerr= a['D18_sterr'].values*2, yerr = np.sqrt(a_gs['D18_sterr']**2+a['D18_sterr']**2).values*2, fmt = 'bs' )
ax.plot(np.linspace(a['D18_stochastic'].min(), a['D18_stochastic'].max()), np.zeros(50), 'g-')
ax.set_xlabel(ur'$\Delta\mathrm{18_{stochastic, adductLine}}$ (\u2030)')
ax.set_ylabel(ur'$\Delta \Delta \mathrm{18_{goalSeek - adductLine}}$ (\u2030)')
plt.savefig('D18_gs_adduct_comparison.pdf')

fig, ax = plt.subplots()
ax.errorbar(a_gs['d13C (vpdb)'] - a['d13C (vpdb)'],a_gs['D18_stochastic'] - a['D18_stochastic'], xerr= np.sqrt(a_gs['d13C_sterr']**2+a['d13C_sterr']**2).values*2, yerr = np.sqrt(a_gs['D18_sterr']**2+a['D18_sterr']**2).values*2, fmt = 'bs' )
ax.set_xlabel(ur'$\Delta\mathrm{^{13}C_{goalSeek - adductLine}}$ (\u2030)')
ax.set_ylabel(ur'$\Delta \Delta \mathrm{18_{goalSeek - adductLine}}$ (\u2030)')
plt.savefig('D18_d13C_gs_adduct_comparison.pdf')

fig, ax = plt.subplots()
ax.errorbar(a['d18'],a_gs['d18'] - a['d18'], xerr= a['d18_sterr'].values*2, yerr = np.sqrt(a_gs['d18_sterr']**2+a['d18_sterr']**2).values*2, fmt = 'bs' )
ax.plot(np.linspace(a['d18'].min(), a['d18'].max()), np.zeros(50), 'g-')
ax.set_xlabel(ur'$\delta\mathrm{18_{stochastic, adductLine}}$ (\u2030)')
ax.set_ylabel(ur'$\Delta\mathrm{18_{goalSeek - adductLine}}$ (\u2030)')
plt.savefig('d18_gs_adduct_comparison_little.pdf')

fig, ax = plt.subplots()
(a_gs['D18_stochastic'] - a['D18_stochastic']).hist(bins = 20)

plt.savefig('D18_histogram.pdf')
