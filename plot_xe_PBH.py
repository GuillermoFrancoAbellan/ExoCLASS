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
                   'Omega_b': 0.0508,
                   'Omega_cdm':0.2645,
                   'H0':67.36,
                   'A_s':2.0989e-9,
                   'n_s':0.9649,
                   'tau_reio':0.0541,
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


pbh_fraction = 0.1
M_mono = 32
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


# EXTENDED PBH, Mcut =100 Msun

if True:

    M = Class()
    M.set(common_settings)
    M.set({
            'recombination': 'recfast',
            'PBH_fraction': pbh_fraction,
            'PBH_accretion_recipe': 'disk_accretion',
            'log10_Mcut_PBH': 2,
            'has_extended_PBH_MassFunc':  'yes',
            'energy_deposition_function': 'DarkAges',
            'DarkAges_mode': 'built_in'
            })
    M.compute()

    thermo = M.get_thermodynamics()
    z_ext1 = thermo['z']
    xe_ext1 = thermo['x_e']
    clM = M.lensed_cl(2500)
    ll_ext1 = clM['ell'][2:]
    clTT_ext1 = clM['tt'][2:]
    M.struct_cleanup()
    M.empty()    
    fxe_ext1 = interp1d(z_ext1,xe_ext1)
    fclTT_ext1 = interp1d(ll_ext1,clTT_ext1)


# EXTENDED PBH, Mcut = 10000 Msun

    M = Class()
    M.set(common_settings)
    M.set({
            'recombination': 'recfast',
            'PBH_fraction': pbh_fraction,
            'PBH_accretion_recipe': 'disk_accretion',
            'log10_Mcut_PBH': 4.5,
            'has_extended_PBH_MassFunc':  'yes',
            'energy_deposition_function': 'DarkAges',
            'DarkAges_mode': 'built_in'
            })
    M.compute()

    thermo = M.get_thermodynamics()
    z_ext2 = thermo['z']
    xe_ext2 = thermo['x_e']
    clM = M.lensed_cl(2500)
    ll_ext2 = clM['ell'][2:]
    clTT_ext2 = clM['tt'][2:]
    M.struct_cleanup()
    M.empty()    
    fxe_ext2 = interp1d(z_ext2,xe_ext2)
    fclTT_ext2 = interp1d(ll_ext2,clTT_ext2)


# MONOCHROMATIC PBH, Mpeak =1 Msun

    M = Class()
    M.set(common_settings)
    M.set({
    'recombination': 'recfast',
    'PBH_fraction': pbh_fraction,
    'PBH_accretion_recipe': 'disk_accretion',
    'PBH_accreting_mass': 1.0,
    'energy_deposition_function': 'DarkAges',
    'DarkAges_mode': 'built_in'
    })
    M.compute() 
    thermo = M.get_thermodynamics()
    z_mon1 = thermo['z']
    xe_mon1 = thermo['x_e']  
    clM = M.lensed_cl(2500)
    ll_mon1 = clM['ell'][2:]
    clTT_mon1 = clM['tt'][2:]  
    M.struct_cleanup()
    M.empty()
    fxe_mon1 = interp1d(z_mon1,xe_mon1)
    fclTT_mon1 = interp1d(ll_mon1,clTT_mon1)

# MONOCHROMATIC PBH, Mpeak = 15 Msun
    
    M = Class()
    M.set(common_settings)
    M.set({
    'recombination': 'recfast',
    'PBH_fraction': pbh_fraction,
    'PBH_accretion_recipe': 'disk_accretion',
    'PBH_accreting_mass': M_mono,
    'energy_deposition_function': 'DarkAges',
    'DarkAges_mode': 'built_in'
    })
    M.compute()
    thermo = M.get_thermodynamics()
    z_mon2 = thermo['z']
    xe_mon2 = thermo['x_e']
    clM = M.lensed_cl(2500)
    ll_mon2 = clM['ell'][2:]
    clTT_mon2 = clM['tt'][2:]
    M.struct_cleanup()
    M.empty()
    fxe_mon2 = interp1d(z_mon2,xe_mon2)
    fclTT_mon2 = interp1d(ll_mon2,clTT_mon2)


#%%
    
gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
ax_1 = plt.subplot(gs[0])
ax_2 = plt.subplot(gs[1],sharex = ax_1)
plt.subplots_adjust(hspace=0)


ax_2.set_xlabel(r"redshift $z$",fontsize=19)
ax_1.set_ylabel(r"free electron fraction $x_e(z)$",fontsize=19)
ax_2.set_ylabel(r"$\frac{x_e(\mathrm{PBH})}{x_e(\Lambda\mathrm{CDM})}-1$",fontsize=19)


ax_1.loglog(z_lcdm,xe_lcdm,color = 'blue', linestyle='solid', label=r'Standard')
ax_1.loglog(z_mon1,xe_mon1,color = 'green', linestyle='dashed', label=r'Monochromatic, $M_{\mathrm{peak}} =1 \ M_{\odot}$')
ax_1.loglog(z_mon2,xe_mon2,color = 'green', linestyle='dotted', label=r'Monochromatic, $M_{\mathrm{peak}} =%.0f \ M_{\odot}$'%M_mono)
ax_1.loglog(z_ext1,xe_ext1,color = 'red', linestyle='dashed', label=r'Extended, $M_{\mathrm{cut}} =10^2 \ M_{\odot}$')
ax_1.loglog(z_ext2,xe_ext2,color = 'red', linestyle='dotted', label=r'Extended, $M_{\mathrm{cut}} =10^{4.5} \ M_{\odot}$')

ax_2.semilogx(z_lcdm,fxe_mon1(z_lcdm)/fxe_lcdm(z_lcdm)-1.,color = 'green', linestyle='dashed', label=r'Monochromatic, $M_{\mathrm{peak}} =1 \ M_{\odot}$')
ax_2.semilogx(z_lcdm,fxe_mon2(z_lcdm)/fxe_lcdm(z_lcdm)-1.,color = 'green', linestyle='dotted', label=r'Monochromatic, $M_{\mathrm{peak}} =%.0f \ M_{\odot}$'%M_mono)
ax_2.semilogx(z_lcdm,fxe_ext1(z_lcdm)/fxe_lcdm(z_lcdm)-1.,color = 'red', linestyle='dashed', label=r'Extended, $M_{\mathrm{cut}} =10^2 \ M_{\odot}$')
ax_2.semilogx(z_lcdm,fxe_ext2(z_lcdm)/fxe_lcdm(z_lcdm)-1,color = 'red', linestyle='dotted', label=r'Extended, $M_{\mathrm{cut}} =10^{4.5} \ M_{\odot}$')


ax_1.text(1200,4.0e-3,r'$f_{\mathrm{PBH}}=%.1f$'%pbh_fraction)
ax_1.text(1200,2.0e-3,r'disk accretion')
ax_1.text(1200,1.1e-3,r'no DM halos')

ax_1.legend(frameon=False,fontsize =16,loc='upper left',borderaxespad=0.)

ax_2.set_ylim([-3,20])
ax_1.set_ylim([0.00015,2])
ax_1.set_xlim([100,2000])
ax_2.set_xlim([100,2000])

ax_2.set_xticks([100,300, 500, 1000, 3000])
ax_2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())


ax_1.tick_params(axis="y", labelsize=16)
ax_2.tick_params(axis="y", labelsize=16)
ax_2.tick_params(axis="x", labelsize=16)
ax_1.tick_params(axis="x", labelsize=1)

plt.show()

#%%

gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
ax_1 = plt.subplot(gs[0])
ax_2 = plt.subplot(gs[1],sharex = ax_1)
plt.subplots_adjust(hspace=0)


ax_2.set_xlabel(r"multipole $\ell$",fontsize=19)
ax_1.set_ylabel(r"$\ell (\ell+1) C_{\ell}^{TT} / 2 \pi $",fontsize=19)
ax_2.set_ylabel(r"$\frac{C_{\ell}(\mathrm{PBH})}{C_{\ell}(\Lambda\mathrm{CDM})}-1$",fontsize=19)


ax_1.loglog(ll_lcdm,clTT_lcdm*ll_lcdm*(ll_lcdm+1.)/2./np.pi,color = 'blue', linestyle='solid', label=r'Standard')
ax_1.loglog(ll_mon1,clTT_mon1*ll_mon1*(ll_mon1+1.)/2./np.pi,color = 'green', linestyle='dashed', label=r'Monochromatic, $M_{\mathrm{peak}} =1 \ M_{\odot}$')
ax_1.loglog(ll_mon2,clTT_mon2*ll_mon2*(ll_mon2+1.)/2./np.pi,color = 'green', linestyle='dotted', label=r'Monochromatic, $M_{\mathrm{peak}} =%.0f \ M_{\odot}$'%M_mono)
ax_1.loglog(ll_ext1,clTT_ext1*ll_ext1*(ll_ext1+1.)/2./np.pi,color = 'red', linestyle='dashed', label=r'Extended, $M_{\mathrm{cut}} =10^2 \ M_{\odot}$')
ax_1.loglog(ll_ext2,clTT_ext2*ll_ext2*(ll_ext2+1.)/2./np.pi,color = 'red', linestyle='dotted', label=r'Extended, $M_{\mathrm{cut}} =10^{4.5} \ M_{\odot}$')

ax_2.semilogx(ll_lcdm,fclTT_mon1(ll_lcdm)/fclTT_lcdm(ll_lcdm)-1.,color = 'green', linestyle='dashed', label=r'Monochromatic, $M_{\mathrm{peak}} =1 \ M_{\odot}$')
ax_2.semilogx(ll_lcdm,fclTT_mon2(ll_lcdm)/fclTT_lcdm(ll_lcdm)-1,color = 'green', linestyle='dotted', label=r'Monochromatic, $M_{\mathrm{peak}} =%.0f \ M_{\odot}$'%M_mono)
ax_2.semilogx(ll_lcdm,fclTT_ext1(ll_lcdm)/fclTT_lcdm(ll_lcdm)-1,color = 'red', linestyle='dashed', label=r'Extended, $M_{\mathrm{cut}} =10^2 \ M_{\odot}$')
ax_2.semilogx(ll_lcdm,fclTT_ext2(ll_lcdm)/fclTT_lcdm(ll_lcdm)-1,color = 'red', linestyle='dotted', label=r'Extended, $M_{\mathrm{cut}} =10^{4.5} \ M_{\odot}$')


ax_1.legend(frameon=False,fontsize =16,loc='upper left',borderaxespad=0.)

ax_1.text(100,5.0e-11,r'$f_{\mathrm{PBH}}=%.1f$'%pbh_fraction)
ax_1.text(100,4.0e-11,r'disk accretion')
ax_1.text(100,3.2e-11,r'no DM halos')

def binned_cosmic_variance (result,l_ini,width):
    result = 0
    Clb = 0
    for i in range(0,int(width)):
        result += 2/(2*(l_ini+float(i))+1)*(l_ini+float(i))*(l_ini+float(i)+1)*(l_ini+float(i))*(l_ini+float(i)+1)*fclTT_lcdm(l_ini+i)*fclTT_lcdm(l_ini+i)
        Clb += (l_ini+float(i))*(l_ini+float(i)+1)*fclTT_lcdm(l_ini+i)
    return np.sqrt(result)/Clb

l_min = 2.;
l_max = 2400;
n_step = 100.;
j=0.
step = l_min
width= 25

while step < l_max:
        result = 0.0
        if step < 29:
            width = 1
            step = l_min+j*width
            j+=1
            if step == 29:
                j = 0
                l_min = 30
        else:
            width = 30
            step = l_min+j*width
            j+=1

        ax_2.add_patch(patches.Rectangle((int(step), -1*binned_cosmic_variance(result,int(step),width)), width, 2*binned_cosmic_variance(result,int(step),width),color='red',alpha=0.1))




ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])
ax_1.set_ylim([2e-11,1e-9])
ax_2.set_ylim([-0.35,0.35])

ax_1.tick_params(axis="y", labelsize=16)
ax_2.tick_params(axis="y", labelsize=16)
ax_2.tick_params(axis="x", labelsize=16)
ax_1.tick_params(axis="x", labelsize=8)

plt.show()
