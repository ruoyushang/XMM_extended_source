import os
import xspec
import pandas as pd
import matplotlib.pylab as plt


file_dir = '../output_plots/plot_Cas_A/'

data_file = "%s/table_spectrum_0400210101.fits"%(file_dir)
rmf_file = "%s/table_rmf_0400210101.fits"%(file_dir)
arf_file = "%s/table_arf_0400210101.fits"%(file_dir)

data = xspec.Spectrum(data_file,respFile=rmf_file,arfFile=arf_file)

model = xspec.Model("bbody+gaus")
ncomp = len(model.componentNames)
for icomp in model.componentNames:
    print (icomp,eval(f'model.{icomp}.parameterNames'))

line_FeXXV = 6.6
model.gaussian.LineE = line_FeXXV
model.gaussian.LineE.frozen = True

startE = 4.5 # keV
endE = 9.0 # keV
data.notice("all")
data.ignore(f"bad")
data.ignore(f"**-{startE} {endE}-**")

model.show()
xspec.Fit.nIterations = 100
xspec.Fit.query = 'yes'
xspec.Fit.perform()
xspec.Fit.show()
cstat = xspec.Fit.statistic
dof = xspec.Fit.dof
chi2r = xspec.Fit.testStatistic/dof
print (f'chi2r = {chi2r}')

# The idea on how to plot the model components is from the XSPEC facebook group, idea from Andy Beardmore. Saving in QDP file.
# Following this example: https://gist.github.com/ivvv/716799056b4aadc87e41097472edbd20
save_data = '%s/specfit.qdp'%(file_dir)
if (os.path.isfile(save_data)):
    os.remove(save_data)
xspec.Plot.device = '/null'
xspec.Plot.xAxis = "keV"
xspec.Plot.add = True
xspec.Plot.addCommand(f'wd {save_data}')
xspec.Plot("ld")

names = ['e','de','rate','rate_err','total']
for j in range(ncomp):
    names.append(f'model{j}')
#print (names)
df = pd.read_table(save_data,skiprows=3,names=names, delimiter=' ')
print (f'df = {df}')

fig, ax = plt.subplots(figsize=(10,6))
ax.errorbar(df.e, df.rate, xerr=df.de,yerr=df.rate_err,fmt='o-',label='data')
ax.plot(df.e, df.total, color='red',label='Total model',linewidth=3)
for j in range(ncomp):
    ax.plot(df.e, df[f'model{j}'],label=f'{model.componentNames[j]}')
ax.set_xlabel('Energy (keV)')
ax.set_ylabel(r'counts/s/keV')
ax.set_xscale("linear")
ax.set_yscale("linear")
#ax.set_ylim((0.001,0.1))
ax.grid()
ax.legend()
fig.savefig('%s/fit_result.png'%(file_dir),bbox_inches='tight')

print ('Done.')
