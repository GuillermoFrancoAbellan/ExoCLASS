import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.lines as mlines


rc('text', usetex=True)

font = {'size': 16, 'family': 'STIXGeneral'}
axislabelfontsize='large'
matplotlib.rc('font', **font)


##################### READ THE DATA ############################################

if True:

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em6_k1e3_Delos_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    zz = dat[:,0]
    fBz_Delos = interp1d(dat[:,0], dat[:,9])


    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em6_k1e3_GG_spike_zFst_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_GG_spike_zFst = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em6_k1e3_GG_spike_zFavg_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_GG_spike_zFavg = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em6_k1e3_GG_all_zFst_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_GG_all_zFst = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH_A1em6_k1e3_GG_all_zFavg_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fBz_GG_all_zFavg = interp1d(dat[:,0], dat[:,9])




##################### BOOST FACTOR ############################################
fig, ax = plt.subplots()

plt.xlabel(r"Redshift $z$",fontsize=19)
plt.ylabel(r"Boost factor $1+\mathcal{B}_0(z)$",fontsize=19)

plt.loglog(zz,fBz_Delos(zz),color = 'black', linestyle='solid', label=r'BBKS')
plt.loglog(zz,fBz_GG_spike_zFst(zz),color = 'red', linestyle='dashed', label=r'EST, only spike, $Z_f(M_{\star}^{-})$')
plt.loglog(zz,fBz_GG_all_zFst(zz),color = 'blue', linestyle='dashed', label=r'EST, spike+smooth, $Z_f(M_{\star}^{-})$')
plt.loglog(zz,fBz_GG_all_zFavg(zz),color = 'blue', linestyle='dotted', label=r'EST, spike+smooth, $\langle Z_f(M_{\star})\rangle$')
plt.loglog(zz,fBz_GG_spike_zFavg(zz),color = 'red', linestyle='dotted', label=r'EST, only spike, $\langle Z_f(M_{\star})\rangle$')

plt.fill_between(zz, fBz_GG_spike_zFavg(zz), fBz_GG_spike_zFst(zz), color='red',alpha=0.5)
plt.fill_between(zz, fBz_GG_all_zFavg(zz), fBz_GG_all_zFst(zz), color='blue',alpha=0.3)


legend1 = plt.legend(frameon=True,fontsize =16,loc='upper right',borderaxespad=0.6)
legend1.get_frame().set_linewidth(2)

ax.add_artist(legend1)

line1 = mlines.Line2D([], [], color='white', linestyle='', label = r'$\mathcal{A}_{\star}= 10^{-6}, \ k_{\star} = 10^3 \ \mathrm{Mpc}^{-1}$')

legend2 = plt.legend(handles = [line1], loc='lower left', fontsize=16, frameon=True,handletextpad=-0.01, handlelength=0 )
legend2.get_frame().set_linewidth(2)

plt.ylim([0.5,1e12])
plt.xlim([1e-1,2e3])

plt.grid(linewidth=0.3)

plt.tick_params(axis="y", labelsize=18,width=1.5,which = 'major')
plt.tick_params(axis="x", labelsize=18,width=1.5,which = 'major')

plt.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')
plt.tick_params(axis="x", labelsize=18,width=1.5,which = 'minor')

plt.setp(ax.spines.values(), linewidth=1.5)

plt.show()
