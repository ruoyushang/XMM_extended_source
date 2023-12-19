import os
import xspec
import pandas as pd
import matplotlib.pylab as plt
import numpy as np


file_dir = '../output_plots/plot_Cas_A/'
#runID = '0400210101'
runID = '0412180101'

#region = 'r0'
region = 'r1'
#region = 'r2'
#region = 'r3'
#region = 'r4'

# running grppha in python: https://github.com/tenoto/fspec/blob/master/fgrppha.py
inpha =  "%s/table_spectrum_data_%s_%s.fits"%(file_dir,runID,region)
outpha = "%s/table_spectrum_data_%s_%s_grppha.fits"%(file_dir,runID,region)
cmd  = 'rm %s\n' % outpha
cmd += 'grppha<<EOF\n'
cmd += '%s\n' % inpha
cmd += '%s\n' % outpha
cmd += "group min 200\n"
cmd += "exit\n"
cmd += "EOF\n"
os.system(cmd)

data_file = "%s/table_spectrum_data_%s_%s_grppha.fits"%(file_dir,runID,region)
bkg_file = "%s/table_spectrum_bkgd_%s_%s.fits"%(file_dir,runID,region)
rmf_file = "%s/table_rmf_%s_%s.fits"%(file_dir,runID,region)
arf_file = "%s/table_arf_%s_%s.fits"%(file_dir,runID,region)
print (f'data_file = {data_file}')
#exit()

data = xspec.Spectrum(data_file,backFile=bkg_file,respFile=rmf_file,arfFile=arf_file)

#model = xspec.Model("bbody+gaus+gaus+gaus+gaus+gaus+gaus+gaus+gaus")
#model = xspec.Model("vapec")
model = xspec.Model("vnei")
ncomp = len(model.componentNames)
print (f'===================================================================================================')
for icomp in model.componentNames:
    print (icomp,eval(f'model.{icomp}.parameterNames'))
#exit()

#model.vapec.kT.frozen = False
#model.vapec.Redshift.frozen = False
##model.vapec.He.frozen = False
##model.vapec.C.frozen = False
##model.vapec.N.frozen = False
##model.vapec.O.frozen = False
##model.vapec.Ne.frozen = False
##model.vapec.Mg.frozen = False
##model.vapec.Al.frozen = False
##model.vapec.Si.frozen = False
##model.vapec.S.frozen = False
##model.vapec.Ar.frozen = False
#model.vapec.Ca.frozen = False
#model.vapec.Fe.frozen = False
##model.vapec.Ni.frozen = False

model.vnei.kT.frozen = False
model.vnei.Redshift.frozen = False
#model.vnei.H.frozen = False
#model.vnei.He.frozen = False
#model.vnei.C.frozen = False
#model.vnei.N.frozen = False
model.vnei.O.frozen = False
model.vnei.Ne.frozen = False
model.vnei.Mg.frozen = False
model.vnei.Si.frozen = False
model.vnei.S.frozen = False
model.vnei.Ar.frozen = False
model.vnei.Ca.frozen = False
model.vnei.Fe.frozen = False
#model.vnei.Ni.frozen = False

#model.nei.kT.frozen = False
#model.nei.Abundanc.frozen = False

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


#startE = 1.0 # keV
startE = 1.7 # keV
endE = 8.0 # keV
data.notice("all")
data.ignore(f"bad")
data.ignore(f"**-{startE} {endE}-**")

line_name = []
line_energy = []
line_name += ['Fe XXV']
line_energy += [6.6]
line_name += ['Ca XIX']
line_energy += [3.902]
#line_name += ['Ar XVIII']
#line_energy += [3.323]
line_name += ['Ar XVII']
line_energy += [3.140]
#line_name += ['S XVI']
#line_energy += [3.107]
#line_name += ['S XV']
#line_energy += [2.884]
#line_name += ['S XVI']
#line_energy += [2.620]
#line_name += ['Si XIV']
#line_energy += [2.506]
line_name += ['S XV']
line_energy += [2.449]
#line_name += ['Si XIV']
#line_energy += [2.376]
#line_name += ['Si XIII']
#line_energy += [2.294]
#line_name += ['Si XIII']
#line_energy += [2.183]
#line_name += ['Si XIV']
#line_energy += [2.0]
line_name += ['Si XIII']
line_energy += [1.865]
#line_name += ['Mg XII']
#line_energy += [1.473]
#line_name += ['Fe XXI']
#line_energy += [1.314]
#line_name += ['Fe XXIII']
#line_energy += [1.129]
#line_name += ['Fe XXII']
#line_energy += [1.053]
#line_name += ['Fe XXIII']
#line_energy += [1.020]
#line_name += ['Fe XXI']
#line_energy += [1.000]


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
#ax.plot(df.e, df.total, color='red',label='vapec model')
ax.plot(df.e, df.total, color='red',label='vnei model')
#for j in range(ncomp):
#    ax.plot(df.e, df[f'model{j}'],label=f'{model.componentNames[j]}')
max_height = np.max(df.total)
for line in range(0,len(line_name)):
    line_idx = find_closest_index(df.e, line_energy[line])
    text_height = df.total[line_idx] * 0.1
    ax.annotate(line_name[line], xy=(line_energy[line], 0.5*df.total[line_idx]), xytext=(line_energy[line], text_height), arrowprops=arrowprops)
textstr = '\n'.join((
    #r'kT = $%.2f \pm %.2f$ %s' % (model.vapec.kT.values[0],model.vapec.kT.sigma,model.vapec.kT.unit),
    #r'Ca = $%.2f \pm %.2f$ %s' % (model.vapec.Ca.values[0],model.vapec.Ca.sigma,model.vapec.Ca.unit),
    #r'Fe = $%.2f \pm %.2f$ %s' % (model.vapec.Fe.values[0],model.vapec.Fe.sigma,model.vapec.Fe.unit),
    r'kT = $%.2f \pm %.2f$ %s' % (model.vnei.kT.values[0],model.vnei.kT.sigma,model.vnei.kT.unit),
    r'Si = $%.2f \pm %.2f$ %s' % (model.vnei.Si.values[0],model.vnei.Si.sigma,model.vnei.Si.unit),
    r'S  = $%.2f \pm %.2f$ %s' % (model.vnei.S.values[0],model.vnei.S.sigma,model.vnei.S.unit),
    r'Ar = $%.2f \pm %.2f$ %s' % (model.vnei.Ar.values[0],model.vnei.Ar.sigma,model.vnei.Ar.unit),
    r'Ca = $%.2f \pm %.2f$ %s' % (model.vnei.Ca.values[0],model.vnei.Ca.sigma,model.vnei.Ca.unit),
    r'Fe = $%.2f \pm %.2f$ %s' % (model.vnei.Fe.values[0],model.vnei.Fe.sigma,model.vnei.Fe.unit),
    ))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.70, 0.85, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)
ax.set_xlabel('Energy (keV)')
ax.set_ylabel(r'counts/s/keV')
ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_ylim((0.001,0.1))
ax.grid()
ax.legend()
fig.savefig('%s/fit_result_%s.png'%(file_dir,region),bbox_inches='tight')

print ('Done.')
