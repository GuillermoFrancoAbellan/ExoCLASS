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

lTT,DlTT_mean,DlTT_error_minus,DlTT_error_plus,DlTT_bestfit= np.loadtxt("error_Planck/Planck2018_errorTT.txt",unpack=True)
lEE,DlEE_mean,DlEE_error_minus,DlEE_error_plus,DlEE_bestfit= np.loadtxt("error_Planck/Planck2018_errorEE.txt",unpack=True)
lTE,DlTE_mean,DlTE_error_minus,DlTE_error_plus,DlTE_bestfit= np.loadtxt("error_Planck/Planck2018_errorTE.txt",unpack=True)


m_DM = 50.0

Log10_A_spike_1 = -7.0
Log10_k_spike_1 = 3.0

Log10_A_spike_2 = -6.0
Log10_k_spike_2 = 3.0

##################### READ THE DATA ############################################

if True:

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/LCDM_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fxe_lcdm = interp1d(dat[:,0], dat[:,2])
    z_lcdm = dat[:,0]

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/LCDM_cl_lensed.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fclTT_lcdm = interp1d(dat[:,0], dat[:,1])
    fclEE_lcdm = interp1d(dat[:,0], dat[:,2])
    ll_lcdm = dat[:,0]

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/smooth_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fxe_smooth = interp1d(dat[:,0], dat[:,2])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/smooth_cl_lensed.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fclTT_smooth = interp1d(dat[:,0], dat[:,1])
    fclEE_smooth = interp1d(dat[:,0], dat[:,2])


    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH1_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fxe_ucmh1 = interp1d(dat[:,0], dat[:,2])
    fBz_ucmh1 = interp1d(dat[:,0], dat[:,9])


    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH1_cl_lensed.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fclTT_ucmh1 = interp1d(dat[:,0], dat[:,1])
    fclEE_ucmh1 = interp1d(dat[:,0], dat[:,2])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH2_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fxe_ucmh2 = interp1d(dat[:,0], dat[:,2])
    fBz_ucmh2 = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH2_cl_lensed.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fclTT_ucmh2 = interp1d(dat[:,0], dat[:,1])
    fclEE_ucmh2 = interp1d(dat[:,0], dat[:,2])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH3_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fxe_ucmh3 = interp1d(dat[:,0], dat[:,2])
    fBz_ucmh3 = interp1d(dat[:,0], dat[:,9])

    files = ['/Users/gfranco/class_repositories/ExoCLASS/output/UCMH3_cl_lensed.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fclTT_ucmh3 = interp1d(dat[:,0], dat[:,1])
    fclEE_ucmh3 = interp1d(dat[:,0], dat[:,2])


##################### BOOST FACTOR ############################################
fig, ax = plt.subplots()

plt.xlabel(r"Redshift $z$",fontsize=19)
plt.ylabel(r"Boost factor $1+\mathcal{B}_0(z)$",fontsize=19)

plt.loglog(z_lcdm,fBz_ucmh3(z_lcdm),color = 'brown', linestyle='dotted', label=r'No spike')
plt.loglog(z_lcdm,fBz_ucmh1(z_lcdm),color = 'blue', linestyle='dashed', label=r'$\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
plt.loglog(z_lcdm,fBz_ucmh2(z_lcdm),color = 'red', linestyle='dashed', label=r'$\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))


legend1 = plt.legend(frameon=True,fontsize =16,loc='lower left',borderaxespad=0.6)
legend1.get_frame().set_linewidth(2)

ax.add_artist(legend1)

plt.ylim([0.5,1e7])
plt.xlim([1e-1,2e3])

plt.grid(linewidth=0.3)

plt.tick_params(axis="y", labelsize=18,width=1.5,which = 'major')
plt.tick_params(axis="x", labelsize=18,width=1.5,which = 'major')

plt.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')
plt.tick_params(axis="x", labelsize=18,width=1.5,which = 'minor')

plt.setp(ax.spines.values(), linewidth=1.5)


plt.show()



##################### FREE ELECTRON FRACTION ############################################
fig, ax = plt.subplots()

plt.xlabel(r"Redshift $z$",fontsize=19)
plt.ylabel(r"Free electron fraction $x_e(z)$",fontsize=19)

plt.loglog(z_lcdm,fxe_lcdm(z_lcdm),color = 'black', linestyle='solid', label=r'Standard')
plt.loglog(z_lcdm,fxe_smooth(z_lcdm),color = 'green', linestyle='solid', label=r'Smooth background')
plt.loglog(z_lcdm,fxe_ucmh3(z_lcdm),color = 'brown', linestyle='dotted', label=r'Halos, no spike')
plt.loglog(z_lcdm,fxe_ucmh1(z_lcdm),color = 'blue', linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \\ \hspace*{4.35em} k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
plt.loglog(z_lcdm,fxe_ucmh2(z_lcdm),color = 'red', linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \\ \hspace*{4.35em} k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))



legend1 = plt.legend(frameon=True,fontsize =16,loc=(0.02,0.4),borderaxespad=0.2)
legend1.get_frame().set_linewidth(2)

#plt.legend(title=r'$\chi \bar{\chi} \ \rightarrow \ b \bar{b}$',loc='upper center')

line1 = mlines.Line2D([], [], color='white', linestyle='', label = r'$\chi \bar{\chi} \ \rightarrow \ b \bar{b}$')
line2 = mlines.Line2D([], [], color='white', linestyle='', label = r'$\sigma_0 c = 3 \times 10^{-26} \ \mathrm{cm}^3/s$')
line3 = mlines.Line2D([], [], color='white', linestyle='', label = r'$ m_{\chi} = %.0f \ \mathrm{GeV}$'%m_DM)


legend2 = plt.legend(handles= [line1,line2,line3], loc='lower left', fontsize=16, frameon=True,handletextpad=-0.01, handlelength=0 )
legend2.get_frame().set_linewidth(2)


ax.add_artist(legend1)
ax.add_artist(legend2)

plt.grid(linewidth=0.3)

plt.ylim([0.00015,2])
plt.xlim([0.1,2000])

plt.tick_params(axis="y", labelsize=18,width=1.5,which = 'major')
plt.tick_params(axis="x", labelsize=18,width=1.5,which = 'major')

plt.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')
plt.tick_params(axis="x", labelsize=18,width=1.5,which = 'minor')


plt.setp(ax.spines.values(), linewidth=1.5)

plt.show()

##################### TT SPECTRUM ############################################

gs = gridspec.GridSpec(2, 1,height_ratios=[2,1])
ax_1 = plt.subplot(gs[0])
ax_2 = plt.subplot(gs[1],sharex = ax_1)
plt.subplots_adjust(hspace=0)


ax_2.set_xlabel(r"Multipole $\ell$",fontsize=19)
ax_1.set_ylabel(r"$\ell (\ell+1) C_{\ell}^{\mathrm{TT}} / 2 \pi $",fontsize=19)
ax_2.set_ylabel(r"$\Delta C^{\mathrm{TT}}_{\ell}/C^{\mathrm{TT}}_{\ell,\Lambda\mathrm{CDM}}$",fontsize=19)

ax_1.loglog(ll_lcdm,fclTT_lcdm(ll_lcdm),color = 'black', linestyle='solid', label=r'Standard')
ax_1.loglog(ll_lcdm,fclTT_smooth(ll_lcdm),color = 'green', linestyle='solid', label=r'Smooth background')
ax_1.loglog(ll_lcdm,fclTT_ucmh3(ll_lcdm),color = 'brown', linestyle='dotted', label=r'Halos, no spike')
ax_1.loglog(ll_lcdm,fclTT_ucmh1(ll_lcdm),color = 'blue',linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
ax_1.loglog(ll_lcdm,fclTT_ucmh2(ll_lcdm),color = 'red', linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))

ax_2.semilogx(ll_lcdm,fclTT_smooth(ll_lcdm)/fclTT_lcdm(ll_lcdm)-1.,color = 'green', linestyle='solid', label=r'Smooth background')
ax_2.semilogx(ll_lcdm,fclTT_ucmh3(ll_lcdm)/fclTT_lcdm(ll_lcdm)-1,color = 'brown', linestyle='dotted', label=r'Halos, no spike')
ax_2.semilogx(ll_lcdm,fclTT_ucmh1(ll_lcdm)/fclTT_lcdm(ll_lcdm)-1,color = 'blue', linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
ax_2.semilogx(ll_lcdm,fclTT_ucmh2(ll_lcdm)/fclTT_lcdm(ll_lcdm)-1,color = 'red', linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))


legend1 = ax_1.legend(frameon=True,fontsize =16,loc='lower left',borderaxespad=0.6)
legend1.get_frame().set_linewidth(2)

line1 = mlines.Line2D([], [], color='white', linestyle='', label = r'$\chi \bar{\chi} \ \rightarrow \ b \bar{b}$')
line2 = mlines.Line2D([], [], color='white', linestyle='', label = r'$\sigma_0 c = 3 \times 10^{-26} \ \mathrm{cm}^3/s$')
line3 = mlines.Line2D([], [], color='white', linestyle='', label = r'$ m_{\chi} = %.0f \ \mathrm{GeV}$'%m_DM)


legend2 = ax_1.legend(handles= [line1,line2,line3], loc='upper left', fontsize=16, frameon=True,handletextpad=-0.01, handlelength=0 )
legend2.get_frame().set_linewidth(2)

ax_1.add_artist(legend1)
ax_1.add_artist(legend2)

l_cosmic_variance = np.linspace(0,48,1000)
l_cosmic_variance_1 = np.linspace(0,30,1000)
l_cosmic_variance_2 = np.linspace(30,48,2)
slope =np.array([0.13,0.0343])

ax_2.fill_between(l_cosmic_variance_1, -0.13,0.13, color='lightgray' )
ax_2.fill_between(l_cosmic_variance_2, -slope, slope, color='lightgray' )
ax_2.fill_between(lTT, -(DlTT_error_plus)/DlTT_mean, +(DlTT_error_plus)/DlTT_mean, color='lightgray')


ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])
ax_1.set_ylim([8e-12,1e-9])
ax_2.set_ylim([-0.12,0.01])


ax_2.grid(linewidth=0.3)

ax_1.tick_params(axis="y", labelsize=18,width=1.5,which = 'major')
ax_1.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')
ax_1.tick_params(axis="x", labelsize=8,width=1.5,which = 'major')
ax_1.tick_params(axis="x", labelsize=8,width=1.5,which = 'minor')


ax_2.tick_params(axis="y", labelsize=18,width=1.5,which = 'major')
ax_2.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')
ax_2.tick_params(axis="x", labelsize=18,width=1.5,which = 'major')
ax_2.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')


plt.setp(ax_1.spines.values(), linewidth=1.5)
plt.setp(ax_2.spines.values(), linewidth=1.5)


plt.show()

##################### EE SPECTRUM ############################################


gs = gridspec.GridSpec(2, 1,height_ratios=[2,1])
ax_1 = plt.subplot(gs[0])
ax_2 = plt.subplot(gs[1],sharex = ax_1)
plt.subplots_adjust(hspace=0)


ax_2.set_xlabel(r"Multipole $\ell$",fontsize=19)
ax_1.set_ylabel(r"$\ell (\ell+1) C_{\ell}^{\mathrm{EE}} / 2 \pi $",fontsize=19)
ax_2.set_ylabel(r"$\Delta C^{\mathrm{EE}}_{\ell}/C^{\mathrm{EE}}_{\ell,\Lambda\mathrm{CDM}}$",fontsize=19)

ax_1.loglog(ll_lcdm,fclEE_lcdm(ll_lcdm),color = 'black', linestyle='solid', label=r'Standard')
ax_1.loglog(ll_lcdm,fclEE_smooth(ll_lcdm),color = 'green', linestyle='solid', label=r'Smooth background')
ax_1.loglog(ll_lcdm,fclEE_ucmh3(ll_lcdm),color = 'brown', linestyle='dotted', label=r'Halos, no spike')
ax_1.loglog(ll_lcdm,fclEE_ucmh1(ll_lcdm),color = 'blue',linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
ax_1.loglog(ll_lcdm,fclEE_ucmh2(ll_lcdm),color = 'red', linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))

ax_2.semilogx(ll_lcdm,fclEE_smooth(ll_lcdm)/fclEE_lcdm(ll_lcdm)-1.,color = 'green', linestyle='solid', label=r'Smooth background')
ax_2.semilogx(ll_lcdm,fclEE_ucmh3(ll_lcdm)/fclEE_lcdm(ll_lcdm)-1,color = 'brown', linestyle='dotted', label=r'Halos, no spike')
ax_2.semilogx(ll_lcdm,fclEE_ucmh1(ll_lcdm)/fclEE_lcdm(ll_lcdm)-1,color = 'blue', linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
ax_2.semilogx(ll_lcdm,fclEE_ucmh2(ll_lcdm)/fclEE_lcdm(ll_lcdm)-1,color = 'red', linestyle='dashed', label=r'Halos, $\mathcal{A}_{\star} = 10^{%.0f}, \ k_{\star} = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))



legend1 = ax_1.legend(frameon=True,fontsize =16,loc='upper left',borderaxespad=0.6)
legend1.get_frame().set_linewidth(2)

line1 = mlines.Line2D([], [], color='white', linestyle='', label = r'$\chi \bar{\chi} \ \rightarrow \ b \bar{b}$')
line2 = mlines.Line2D([], [], color='white', linestyle='', label = r'$\sigma_0 c = 3 \times 10^{-26} \ \mathrm{cm}^3/s$')
line3 = mlines.Line2D([], [], color='white', linestyle='', label = r'$ m_{\chi} = %.0f \ \mathrm{GeV}$'%m_DM)

legend2 = ax_1.legend(handles= [line1,line2,line3], loc='lower right', fontsize=16, frameon=True,handletextpad=-0.01, handlelength=0 )
legend2.get_frame().set_linewidth(2)

ax_1.add_artist(legend1)
ax_1.add_artist(legend2)


l_cosmic_variance = np.linspace(0,48,1000)
l_cosmic_variance_1 = np.linspace(0,30,1000)
l_cosmic_variance_2 = np.linspace(30,48,2)
slope =np.array([0.13,0.0343])

ax_2.fill_between(l_cosmic_variance, -0.26,0.26, color='lightgray' )

ax_2.fill_between(lEE, -(DlEE_error_plus)/DlEE_mean, +(DlEE_error_plus)/DlEE_mean, color='lightgray')


ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])
ax_1.set_ylim([2e-16,6e-11])
ax_2.set_ylim([-0.25,0.25])

ax_2.grid(linewidth=0.3)

ax_1.tick_params(axis="y", labelsize=18,width=1.5,which = 'major')
ax_1.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')
ax_1.tick_params(axis="x", labelsize=8,width=1.5,which = 'major')
ax_1.tick_params(axis="x", labelsize=8,width=1.5,which = 'minor')


ax_2.tick_params(axis="y", labelsize=18,width=1.5,which = 'major')
ax_2.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')
ax_2.tick_params(axis="x", labelsize=18,width=1.5,which = 'major')
ax_2.tick_params(axis="y", labelsize=18,width=1.5,which = 'minor')


plt.setp(ax_1.spines.values(), linewidth=1.5)
plt.setp(ax_2.spines.values(), linewidth=1.5)


plt.show()
