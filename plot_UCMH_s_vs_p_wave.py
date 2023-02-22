import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.lines as mlines
import itertools
import seaborn as sns

rc('text', usetex=True)

font = {'size': 16, 'family': 'STIXGeneral'}
axislabelfontsize='large'
matplotlib.rc('font', **font)


##################### READ THE DATA ############################################

if True:

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em5_k1e0_s_wave_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    zz = dat[:,0]
    fBz_A1em5_k1e0_s_wave = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em5_k1e1_s_wave_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_A1em5_k1e1_s_wave = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em5_k1e2_s_wave_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_A1em5_k1e2_s_wave = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em5_k1e3_s_wave_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_A1em5_k1e3_s_wave = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em5_k1e0_p_wave_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    zz = dat[:,0]
    fBz_A1em5_k1e0_p_wave = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em5_k1e1_p_wave_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_A1em5_k1e1_p_wave = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em5_k1e2_p_wave_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_A1em5_k1e2_p_wave = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em5_k1e3_p_wave_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_A1em5_k1e3_p_wave = interp1d(dat[:,0], dat[:,9])


##################### BOOST FACTOR ############################################
fig, ax = plt.subplots()

plt.xlabel(r"Redshift $z$",fontsize=19)
plt.ylabel(r"Boost factor $\mathcal{R}_\ell(z)$",fontsize=19)

pal = itertools.cycle(sns.color_palette("viridis"))

plt.loglog(zz,fBz_A1em5_k1e0_s_wave(zz),color = next(pal), linestyle='solid', label=r'$k_{\star} = 1 \ \mathrm{Mpc}^{-1}$')
next(pal)
plt.loglog(zz,fBz_A1em5_k1e1_s_wave(zz),color = next(pal), linestyle='solid', label=r'$k_{\star} = 10 \ \mathrm{Mpc}^{-1}$')
next(pal)
plt.loglog(zz,fBz_A1em5_k1e2_s_wave(zz),color = next(pal), linestyle='solid', label=r'$k_{\star} = 10^2 \ \mathrm{Mpc}^{-1}$')
plt.loglog(zz,fBz_A1em5_k1e3_s_wave(zz),color = next(pal), linestyle='solid', label=r'$k_{\star} = 10^3 \ \mathrm{Mpc}^{-1}$')

pal = itertools.cycle(sns.color_palette("viridis"))

plt.loglog(zz,fBz_A1em5_k1e0_p_wave(zz),color = next(pal), linestyle='dashdot')
next(pal)
plt.loglog(zz,fBz_A1em5_k1e1_p_wave(zz),color = next(pal), linestyle='dashdot')
next(pal)
plt.loglog(zz,fBz_A1em5_k1e2_p_wave(zz),color = next(pal), linestyle='dashdot')
plt.loglog(zz,fBz_A1em5_k1e3_p_wave(zz),color = next(pal), linestyle='dashdot')


legend1 = plt.legend(frameon=True,fontsize =16,loc='upper right',borderaxespad=0.6)
legend1.get_frame().set_linewidth(2)

ax.add_artist(legend1)

line1 = mlines.Line2D([], [], color='black', linestyle='solid', label   = r's-wave, $ \ \mathcal{A}_{\star}= 10^{-5}$')
line2 = mlines.Line2D([], [], color='black', linestyle='dashdot', label = r'p-wave, $ \ \mathcal{A}_{\star}= 10^{-5}$')

legend2 = plt.legend(handles = [line1, line2], loc='center left', fontsize=16, frameon=True)
legend2.get_frame().set_linewidth(2)

plt.ylim([1e-5,1e8])
plt.xlim([1e-1,5e3])

plt.grid(linewidth=0.3)

plt.tick_params(axis="y", labelsize=18,width=1.5,which = 'major')
plt.tick_params(axis="x", labelsize=18,width=1.5,which = 'major')

plt.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')
plt.tick_params(axis="x", labelsize=18,width=1.5,which = 'minor')

plt.setp(ax.spines.values(), linewidth=1.5)

plt.show()
