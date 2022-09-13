import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from classy import Class
from scipy.interpolate import interp1d
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

rc('font',**{'family':'serif','serif':['Times']}, **{'weight':'heavy'})
rc('text', usetex=True)


font = {'size': 16, 'family': 'STIXGeneral'}
axislabelfontsize='large'
matplotlib.rc('font', **font)

common_settings = {'write background':'yes',
                   'output':'tCl,pCl,lCl',
                   'lensing':'yes',
                   'omega_b': 0.02233,
                   'omega_cdm':0.1198,
                   'H0':67.36,
                   'ln10^{10}A_s':3.043,
                   'n_s':0.9652,
                   'tau_reio':0.0540,
                   'N_ncdm':1,
                   'N_ur': 2.0328,
                   'm_ncdm': 0.06,
                   'T_ncdm': 0.71611,
                   'input_verbose':1,
                   'background_verbose':1,
                   'thermodynamics_verbose':1,
                   'perturbations_verbose':1,
                   'transfer_verbose':1,
                   'primordial_verbose':1,
                   'spectra_verbose': 1,
                   'nonlinear_verbose':1,
                   'lensing_verbose':1,
                   'output_verbose':1
                   }

lTT,DlTT_mean,DlTT_error_minus,DlTT_error_plus,DlTT_bestfit= np.loadtxt("error_Planck/Planck2018_errorTT.txt",unpack=True)
lEE,DlEE_mean,DlEE_error_minus,DlEE_error_plus,DlEE_bestfit= np.loadtxt("error_Planck/Planck2018_errorEE.txt",unpack=True)
lTE,DlTE_mean,DlTE_error_minus,DlTE_error_plus,DlTE_bestfit= np.loadtxt("error_Planck/Planck2018_errorTE.txt",unpack=True)


#%%

# LCDM 
M = Class()
M.set(common_settings)
M.compute()

thermo = M.get_thermodynamics()
z_lcdm = thermo['z']
xe_lcdm = thermo['x_e']

clM = M.lensed_cl(2500)
ll_lcdm = clM['ell'][2:]
clTT_lcdm = clM['tt'][2:]

M.struct_cleanup()
M.empty()

fxe_lcdm = interp1d(z_lcdm,xe_lcdm)
fclTT_lcdm = interp1d(ll_lcdm,clTT_lcdm)


m_DM = 20.0

Log10_A_spike_1 = -7.0
Log10_k_spike_1 = 3.0
    
Log10_A_spike_2 = -6.0
Log10_k_spike_2 = 3.0


if True:

# ANNIHILATION in SMOOTH BACKGROUND

    M = Class()
    M.set(common_settings)
    M.set({
            'recombination': 'hyrec',
            'DM_mass':m_DM,  #in GeV units
            'annihilation_cross_section':3.0e-26,  #in cm^3/s
            'energy_deposition_function': 'DarkAges',
            'DarkAges_mode': 'built_in',
            'injected_particle_spectra':'bottom',
            'injected_particle_branching_ratio': 1
            })
    M.compute()

    thermo = M.get_thermodynamics()
    z_smooth = thermo['z']
    xe_smooth = thermo['x_e']
    clM = M.lensed_cl(2500)
    ll_smooth = clM['ell'][2:]
    clTT_smooth = clM['tt'][2:]
    
    M.struct_cleanup()
    M.empty()  
    
    fxe_smooth = interp1d(z_smooth,xe_smooth)
    fclTT_smooth = interp1d(ll_smooth,clTT_smooth)

    
if False:    
# ANNIHILATION WITHIN UCMHS, MODEL 1

    M = Class()
    M.set(common_settings)
    M.set({
            'recombination': 'hyrec',
            'DM_mass':m_DM,  #in GeV units
            'annihilation_cross_section':3e-26,  #in cm^3/s
            'has_UCMH_spike':'yes',
            'Log10_A_spike':Log10_A_spike_1,
            'Log10_k_spike':Log10_k_spike_1,
            'Mass_min': 1e-9,
            'energy_deposition_function': 'DarkAges',
            'DarkAges_mode': 'built_in',
            'injected_particle_spectra':'bottom',
            'injected_particle_branching_ratio': 1
            })
    M.compute()

    thermo = M.get_thermodynamics()
    z_ucmh1 = thermo['z']
    xe_ucmh1 = thermo['x_e']
    clM = M.lensed_cl(2500)
    ll_ucmh1 = clM['ell'][2:]
    clTT_ucmh1 = clM['tt'][2:]
    
    M.struct_cleanup()
    M.empty()    
    fxe_ucmh1 = interp1d(z_ucmh1,xe_ucmh1)
    fclTT_ucmh1 = interp1d(ll_ucmh1,clTT_ucmh1)

# ANNIHILATION WITHIN UCMHS, MODEL 2

    M = Class()
    M.set(common_settings)
    M.set({
            'recombination': 'hyrec',
            'DM_mass':m_DM,  #in GeV units
            'annihilation_cross_section':3.0e-26,  #in cm^3/s
            'has_UCMH_spike':'yes',
            'Log10_A_spike':Log10_A_spike_2,
            'Log10_k_spike':Log10_k_spike_2,
            'Mass_min': 1e-9,
            'energy_deposition_function': 'DarkAges',
            'DarkAges_mode': 'built_in',
            'injected_particle_spectra':'bottom',
            'injected_particle_branching_ratio': 1

    })
    M.compute() 
    thermo = M.get_thermodynamics()
    z_ucmh2 = thermo['z']
    xe_ucmh2 = thermo['x_e']  
    clM = M.lensed_cl(2500)
    ll_ucmh2 = clM['ell'][2:]
    clTT_ucmh2 = clM['tt'][2:]  
    
    M.struct_cleanup()
    M.empty()
    fxe_ucmh2 = interp1d(z_ucmh2,xe_ucmh2)
    fclTT_ucmh2 = interp1d(ll_ucmh2,clTT_ucmh2)
    
if True:
    files = ['/Users/gfranco/ExoCLASS/output/UCMH1_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fxe_ucmh1 = interp1d(dat[:,0], dat[:,2])
    fBz_ucmh1 = interp1d(dat[:,0], dat[:,9])

    
    files = ['/Users/gfranco/ExoCLASS/output/UCMH1_cl_lensed.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fclTT_ucmh1 = interp1d(dat[:,0], dat[:,1])
       
    files = ['/Users/gfranco/ExoCLASS/output/UCMH2_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fxe_ucmh2 = interp1d(dat[:,0], dat[:,2])
    fBz_ucmh2 = interp1d(dat[:,0], dat[:,9])
    
    files = ['/Users/gfranco/ExoCLASS/output/UCMH2_cl_lensed.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fclTT_ucmh2 = interp1d(dat[:,0], dat[:,1])
    
    files = ['/Users/gfranco/ExoCLASS/output/UCMH3_thermodynamics.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fxe_ucmh3 = interp1d(dat[:,0], dat[:,2])
    fBz_ucmh3 = interp1d(dat[:,0], dat[:,9])
    
    files = ['/Users/gfranco/ExoCLASS/output/UCMH3_cl_lensed.dat']
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))
    dat = data[0]
    fclTT_ucmh3 = interp1d(dat[:,0], dat[:,1])


#%%
    
#gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
#ax_1 = plt.subplot(gs[0])
#ax_2 = plt.subplot(gs[1],sharex = ax_1)
#plt.subplots_adjust(hspace=0)


plt.xlabel(r"redshift $z$",fontsize=19)
plt.ylabel(r"free electron fraction $x_e(z)$",fontsize=19)
#plt.ylabel(r"$\Delta x_e/x_{e,\Lambda\mathrm{CDM}}$",fontsize=19)


plt.loglog(z_lcdm,fxe_lcdm(z_lcdm),color = 'black', linestyle='solid', label=r'Standard')
plt.loglog(z_lcdm,fxe_smooth(z_lcdm),color = 'green', linestyle='solid', label=r'Smooth background')
plt.loglog(z_lcdm,fxe_ucmh3(z_lcdm),color = 'brown', linestyle='dotted', label=r'Halos, no spike')
plt.loglog(z_lcdm,fxe_ucmh1(z_lcdm),color = 'blue', linestyle='dashed', label=r'Halos, $A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
plt.loglog(z_lcdm,fxe_ucmh2(z_lcdm),color = 'red', linestyle='dashed', label=r'Halos, $A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))

#ax_2.semilogx(z_lcdm,fxe_smooth(z_lcdm)/fxe_lcdm(z_lcdm)-1.,color = 'green', linestyle='solid', label=r'Smooth background')
#ax_2.semilogx(z_lcdm,fxe_ucmh1(z_lcdm)/fxe_lcdm(z_lcdm)-1.,color = 'red', linestyle='dashed', label=r'$A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
#ax_2.semilogx(z_lcdm,fxe_ucmh2(z_lcdm)/fxe_lcdm(z_lcdm)-1.,color = 'purple', linestyle='dashed', label=r'$A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))


plt.text(30,1.0,r'$\chi \bar{\chi} \ \rightarrow \ b \bar{b}$')
plt.text(30,0.53,r'$\langle \sigma v \rangle = 3 \times 10^{-26} \ \mathrm{cm}^3/s$')
plt.text(30,0.3,r'$m_{\mathrm{DM}} = %.0f \ \mathrm{GeV}$'%m_DM)

plt.legend(frameon=False,fontsize =16,loc='lower left',borderaxespad=0.)

#ax_2.ylim([-3,20])
plt.ylim([0.00015,2])
plt.xlim([0.1,2000])

#ax_2.xlim([0.1,2000])

#ax_2.set_xticks([10,100,300, 500, 1000, 3000])
#ax_2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())


plt.tick_params(axis="y", labelsize=16)
#ax_2.tick_params(axis="y", labelsize=16)
#ax_2.tick_params(axis="x", labelsize=16)
plt.tick_params(axis="x", labelsize=16)

plt.show()

#%%

gs = gridspec.GridSpec(2, 1,height_ratios=[2,1])
ax_1 = plt.subplot(gs[0])
ax_2 = plt.subplot(gs[1],sharex = ax_1)
plt.subplots_adjust(hspace=0)


ax_2.set_xlabel(r"multipole $\ell$",fontsize=19)
ax_1.set_ylabel(r"$\ell (\ell+1) C_{\ell}^{\mathrm{TT}} / 2 \pi $",fontsize=19)
ax_2.set_ylabel(r"$\Delta C^{\mathrm{TT}}_{\ell}/C^{\mathrm{TT}}_{\ell,\Lambda\mathrm{CDM}}$",fontsize=19)

ax_1.loglog(ll_lcdm,fclTT_lcdm(ll_lcdm)*(ll_lcdm*(ll_lcdm+1.)/(2.*np.pi)),color = 'black', linestyle='solid', label=r'Standard')
ax_1.loglog(ll_lcdm,fclTT_smooth(ll_lcdm)*(ll_lcdm*(ll_lcdm+1.)/(2.*np.pi)),color = 'green', linestyle='solid', label=r'Smooth background')
ax_1.loglog(ll_lcdm,fclTT_ucmh3(ll_lcdm),color = 'brown', linestyle='dotted', label=r'Halos, no spike')
ax_1.loglog(ll_lcdm,fclTT_ucmh1(ll_lcdm),color = 'blue',linestyle='dashed', label=r'Halos, $A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
ax_1.loglog(ll_lcdm,fclTT_ucmh2(ll_lcdm),color = 'red', linestyle='dashed', label=r'Halos, $A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))

ax_2.semilogx(ll_lcdm,fclTT_smooth(ll_lcdm)/fclTT_lcdm(ll_lcdm)-1.,color = 'green', linestyle='solid', label=r'Smooth background')
ax_2.semilogx(ll_lcdm,fclTT_ucmh3(ll_lcdm)/((ll_lcdm*(ll_lcdm+1.)/(2.*np.pi))*fclTT_lcdm(ll_lcdm))-1,color = 'brown', linestyle='dotted', label=r'Halos, no spike')
ax_2.semilogx(ll_lcdm,fclTT_ucmh1(ll_lcdm)/((ll_lcdm*(ll_lcdm+1.)/(2.*np.pi))*fclTT_lcdm(ll_lcdm))-1,color = 'blue', linestyle='dashed', label=r'Halos, $A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
ax_2.semilogx(ll_lcdm,fclTT_ucmh2(ll_lcdm)/((ll_lcdm*(ll_lcdm+1.)/(2.*np.pi))*fclTT_lcdm(ll_lcdm))-1,color = 'red', linestyle='dashed', label=r'Halos, $A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))


ax_1.legend(frameon=False,fontsize =16,loc='lower left',borderaxespad=0.)

ax_1.text(3.5,7.8e-10,r'$\chi \bar{\chi} \ \rightarrow \ b \bar{b}$')
ax_1.text(3.5,5.5e-10,r'$\langle \sigma v \rangle = 3 \times 10^{-26} \ \mathrm{cm}^3/s$')
ax_1.text(3.5,4.0e-10,r'$m_{\mathrm{DM}} = %.0f \ \mathrm{GeV}$'%m_DM)


l_cosmic_variance_1 = np.linspace(0,30,1000)
l_cosmic_variance_2 = np.linspace(30,48,2)
slope =np.array([0.13,0.0343])

ax_2.fill_between(l_cosmic_variance_1, -0.13,0.13, color='lightgray' )
ax_2.fill_between(l_cosmic_variance_2, -slope, slope, color='lightgray' )
ax_2.fill_between(lTT, -(DlTT_error_plus)/DlTT_mean, +(DlTT_error_plus)/DlTT_mean, color='lightgray')


ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])
ax_1.set_ylim([2e-11,1e-9])
ax_2.set_ylim([-0.12,0.01])

ax_1.tick_params(axis="y", labelsize=16)
ax_2.tick_params(axis="y", labelsize=16)
ax_2.tick_params(axis="x", labelsize=16)
ax_1.tick_params(axis="x", labelsize=8)

plt.show()


#%%

plt.xlabel(r"redshift $z$",fontsize=19)
plt.ylabel(r"Boost factor $B(z)$",fontsize=19)

plt.loglog(z_lcdm,fBz_ucmh3(z_lcdm),color = 'brown', linestyle='dotted', label=r'No spike')
plt.loglog(z_lcdm,fBz_ucmh1(z_lcdm),color = 'blue', linestyle='dashed', label=r'$A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_1,Log10_k_spike_1))
plt.loglog(z_lcdm,fBz_ucmh2(z_lcdm),color = 'red', linestyle='dashed', label=r'$A_0 = 10^{%.0f}, \ k_s = 10^{%.0f} \ \mathrm{Mpc}^{-1}$'%(Log10_A_spike_2,Log10_k_spike_2))


plt.legend(frameon=False,fontsize =16,loc='upper right',borderaxespad=0.)


plt.ylim([0.5,1e8])
plt.xlim([1e-2,1e4])

plt.tick_params(axis="y", labelsize=16)
plt.tick_params(axis="x", labelsize=16)

plt.show()
