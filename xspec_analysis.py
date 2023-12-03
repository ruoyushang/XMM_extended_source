import os
import xspec
import pandas as pd
import matplotlib.pylab as plt
import numpy as np


file_dir = '../output_plots/plot_Cas_A/'
#runID = '0400210101'
runID = '0412180101'

data_file = "%s/table_spectrum_%s.fits"%(file_dir,runID)
rmf_file = "%s/table_rmf_%s.fits"%(file_dir,runID)
arf_file = "%s/table_arf_%s.fits"%(file_dir,runID)
print (f'data_file = {data_file}')
#exit()

data = xspec.Spectrum(data_file,respFile=rmf_file,arfFile=arf_file)

#model = xspec.Model("bbody+gaus+gaus+gaus+gaus+gaus+gaus+gaus+gaus")
model = xspec.Model("vapec")
#model = xspec.Model("vvapec")
ncomp = len(model.componentNames)
print (f'===================================================================================================')
for icomp in model.componentNames:
    print (icomp,eval(f'model.{icomp}.parameterNames'))
#exit()

model.vapec.kT.frozen = False
model.vapec.He.frozen = False
model.vapec.C.frozen = False
model.vapec.N.frozen = False
model.vapec.O.frozen = False
model.vapec.Ne.frozen = False
model.vapec.Mg.frozen = False
model.vapec.Al.frozen = False
model.vapec.Si.frozen = False
model.vapec.S.frozen = False
model.vapec.Ar.frozen = False
model.vapec.Ca.frozen = False
model.vapec.Fe.frozen = False
model.vapec.Ni.frozen = False

#model.vvapec.kT.frozen = False
#model.vvapec.H.frozen = False
#model.vvapec.He.frozen = False
#model.vvapec.Li.frozen = False
#model.vvapec.Be.frozen = False
#model.vvapec.B.frozen = False
#model.vvapec.C.frozen = False
#model.vvapec.N.frozen = False
#model.vvapec.O.frozen = False
#model.vvapec.F.frozen = False
#model.vvapec.Ne.frozen = False
#model.vvapec.Na.frozen = False
#model.vvapec.Mg.frozen = False
#model.vvapec.Al.frozen = False
#model.vvapec.Si.frozen = False
#model.vvapec.P.frozen = False
#model.vvapec.S.frozen = False
#model.vvapec.Cl.frozen = False
#model.vvapec.Ar.frozen = False
#model.vvapec.K.frozen = False
#model.vvapec.Ca.frozen = False
#model.vvapec.Sc.frozen = False
#model.vvapec.Ti.frozen = False
#model.vvapec.V.frozen = False
#model.vvapec.Cr.frozen = False
#model.vvapec.Mn.frozen = False
#model.vvapec.Fe.frozen = False
#model.vvapec.Co.frozen = False
#model.vvapec.Ni.frozen = False
#model.vvapec.Cu.frozen = False
#model.vvapec.Zn.frozen = False

#model.bbody.kT.frozen = False
#line_FeXXV = 6.6
#model.gaussian.LineE = line_FeXXV
#model.gaussian.LineE.frozen = True
#line_CaXIX = 3.84
#model.gaussian_3.LineE = line_CaXIX
#model.gaussian_3.LineE.frozen = True
#line_ArXVII = 3.08
#model.gaussian_4.LineE = line_ArXVII
#model.gaussian_4.LineE.frozen = True
#line_SXV_a = 3.00
#model.gaussian_5.LineE = line_SXV_a
#model.gaussian_5.LineE.frozen = True
#line_SXV_b = 2.86
#model.gaussian_6.LineE = line_SXV_b
#model.gaussian_6.LineE.frozen = True
#line_SXV_c = 2.41
#model.gaussian_7.LineE = line_SXV_c
#model.gaussian_7.LineE.frozen = True
#line_SiXIII = 2.16
#model.gaussian_8.LineE = line_SiXIII
#model.gaussian_8.LineE.frozen = True
#line_SiXIV = 2.00
#model.gaussian_9.LineE = line_SiXIV
#model.gaussian_9.LineE.frozen = True


startE = 1.0 # keV
#startE = 3.0 # keV
endE = 9.0 # keV
data.notice("all")
data.ignore(f"bad")
data.ignore(f"**-{startE} {endE}-**")

line_name = []
line_energy = []
line_name += ['Fe XXV']
line_energy += [6.6]
line_name += ['Ca XIX']
line_energy += [3.902]
line_name += ['Ar XVIII']
line_energy += [3.323]
line_name += ['Ar XVII']
line_energy += [3.140]
line_name += ['S XVI']
line_energy += [3.107]
line_name += ['S XV']
line_energy += [2.884]
line_name += ['S XVI']
line_energy += [2.620]
line_name += ['Si XIV']
line_energy += [2.506]
line_name += ['S XV']
line_energy += [2.449]
line_name += ['Si XIV']
line_energy += [2.376]
line_name += ['Si XIII']
line_energy += [2.294]
line_name += ['Si XIII']
line_energy += [2.183]
line_name += ['Si XIV']
line_energy += [2.0]
line_name += ['Si XIII']
line_energy += [1.865]
line_name += ['Mg XII']
line_energy += [1.473]
line_name += ['Fe XXI']
line_energy += [1.314]
line_name += ['Fe XXIII']
line_energy += [1.129]
line_name += ['Fe XXII']
line_energy += [1.053]
line_name += ['Fe XXIII']
line_energy += [1.020]
line_name += ['Fe XXI']
line_energy += [1.000]


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

def find_closest_index(arr, x):
    # Use the min function with a custom key function
    closest_index = min(range(len(arr)), key=lambda i: abs(arr[i] - x))
    return closest_index

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 10,
        'rotation': 90.,
        }
arrowprops = dict(arrowstyle="->")

fig, ax = plt.subplots(figsize=(10,6))
#ax.errorbar(df.e, df.rate, xerr=df.de,yerr=df.rate_err,label='data')
ax.plot(df.e, df.rate, label='data')
ax.plot(df.e, df.total, color='red',label='VApec model (kT = %0.3f +/- %0.3f %s)'%(model.vapec.kT.values[0],model.vapec.kT.sigma,model.vapec.kT.unit))
#for j in range(ncomp):
#    ax.plot(df.e, df[f'model{j}'],label=f'{model.componentNames[j]}')
max_height = np.max(df.total)
for line in range(0,len(line_name)):
    line_idx = find_closest_index(df.e, line_energy[line])
    text_height = df.total[line_idx] + 0.1*(line % 2 + 1)*max_height
    #ax.text(line_energy[line], text_height, line_name[line], fontdict=font)
    ax.annotate(line_name[line], xy=(line_energy[line], df.total[line_idx]), xytext=(line_energy[line], text_height), arrowprops=arrowprops)
ax.set_xlabel('Energy (keV)')
ax.set_ylabel(r'counts/s/keV')
ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_ylim((0.001,0.1))
ax.grid()
ax.legend()
fig.savefig('%s/fit_result.png'%(file_dir),bbox_inches='tight')

print (f'model.vapec.kT.values = {model.vapec.kT.values}')
print (f'model.vapec.kT.sigma = {model.vapec.kT.sigma}')
print (f'model.vapec.kT.unit = {model.vapec.kT.unit}')
print ('Done.')
