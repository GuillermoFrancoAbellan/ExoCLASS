/** @file thermodynamics.c Documented thermodynamics module
 *
 * Julien Lesgourgues, 6.09.2010
 *
 * Deals with the thermodynamical evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the thermodynamics, i.e. to
 *   integrate the thermodynamical equations, and store all
 *   thermodynamical quantities as a function of redshift inside an
 *   interpolation table. The current version of recombination is
 *   based on RECFAST v1.5. The current version of reionization is
 *   based on exactly the same reionization function as in CAMB, in
 *   order to make allow for comparison. It should be easy to
 *   generalize the module to more complicated reionization histories.
 *
 * - to provide a routine which allow other modules to evaluate any
 *   thermodynamical quantities at a given redshift value (by
 *   interpolating within the interpolation table).
 *
 *
 * The logic is the following:
 *
 * - in a first step, the code assumes that there is no reionization,
 *   and computes the ionization fraction, Thomson scattering rate,
 *   baryon temperature, etc., using RECFAST. The result is stored in
 *   a temporary table 'recombination_table' (within a temporary
 *   structure of type 'recombination') for each redshift in a range 0
 *   < z < z_initial.  The sampling in z space is done with a simple
 *   linear step size.
 * - in a second step, the code adds the reionization history,
 *   starting from a redshift z_reio_start. The ionization fraction at
 *   this redshift is read in the previous recombination table in
 *   order to ensure a perfect matching. The code computes the
 *   ionization fraction, Thomson scattering rate, baryon temperature,
 *   etc., using a given parametrization of the reionization
 *   history. The result is stored in a temporary table
 *   'reionization_table' (within a temporary structure of type
 *   'reionization') for each redshift in the range 0 < z <
 *   z_reio_start. The sampling in z space is found automatically,
 *   given the precision parameter 'reionization_sampling'.
 *
 * - in a third step, the code merges the two tables
 *   'recombination_table' and 'reionization_table' inside the table
 *   'thermodynamics_table', and the temporary structures
 *   'recombination' and 'reionization' are freed. In
 *   'thermodynamics_table', the sampling in z space is the one
 *   defined in the recombination algorithm for z_reio_start < z <
 *   z_initial, and the one defined in the reionization algorithm for
 *   0 < z < z_reio_start.
 *
 * - at this stage, only a few columns in the table
 *   'thermodynamics_table' have been filled. In a fourth step, the
 *   remaining columns are filled, using some numerical
 *   integration/derivation routines from the 'array.c' tools module.
 *
 * - small detail: one of the columns contains the maximum variation
 *   rate of a few relevant thermodynamical quantities. This rate
 *   will be used for defining automatically the sampling step size in
 *   the perturbation module. Hence, the exact value of this rate is
 *   unimportant, but its order of magnitude at a given z defines the
 *   sampling precision of the perturbation module. Hence, it is
 *   harmless to use a smoothing routine in order to make this rate
 *   look nicer, although this will not affect the final result
 *   significantly. The last step in the thermodynamics_init module is
 *   to perform this smoothing.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# thermodynamics_init() at the beginning (but after background_init())
 * -# thermodynamics_at_z() at any later time
 * -# thermodynamics_free() at the end, when no more calls to thermodynamics_at_z() are needed
 */

#include "thermodynamics.h"
#include "integration_routines.h" // GFA
#include "root_finding_routines.h"  // GFA
#ifdef HYREC
#include "hyrec.h"
#endif

#ifdef COSMOREC
#include "CosmoRec.h"
#endif

/**
 * Thermodynamics quantities at given redshift z.
 *
 * Evaluates all thermodynamics quantities at a given value of
 * the redshift by reading the pre-computed table and interpolating.
 *
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure (containing pre-computed table)
 * @param z          Input: redshift
 * @param inter_mode Input: interpolation mode (normal or growing_closeby)
 * @param last_index Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback   Input: vector of background quantities (used only in case z>z_initial for getting ddkappa and dddkappa; in that case, should be already allocated and filled, with format short_info or larger; in other cases, will be ignored)
 * @param pvecthermo Output: vector of thermodynamics quantities (assumed to be already allocated)
 * @return the error status
 */

int thermodynamics_at_z(
                        struct background * pba,
                        struct thermo * pth,
                        double z,
                        short inter_mode,
                        int * last_index,
                        double * pvecback,
                        double * pvecthermo
                        ) {

  /** Summary: */

  /** - define local variables */

  double x0;
  double x,Fx,dFx,ddFx, S, dddmu_tilde;
  double a = pba->a_today/(1+z);
	// double tau_at_z;
	// double *pvecback_new;
	// int last_index_new;

  /* - the fact that z is in the pre-computed range 0 <= z <= z_initial
     will be checked in the interpolation routines below. Before
     trying to interpolate, allow the routine to deal with the case z
     > z_intial: then, all relevant quantities can be extrapolated
     using simple analytic approximations */

  if (z >= pth->z_table[pth->tt_size-1]) {
    /* ionization fraction assumed to remain constant at large z */
    x0= pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_xe];
    pvecthermo[pth->index_th_xe] = x0;

    /* Calculate dkappa/dtau (dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T in units of 1/Mpc) */
    pvecthermo[pth->index_th_dkappa] = (1.+z) * (1.+z) * pth->n_e * x0 * _sigma_ * _Mpc_over_m_;

    /* tau_d scales like (1+z)**2 */
    pvecthermo[pth->index_th_tau_d] = pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_d]*pow((1+z)/(1.+pth->z_table[pth->tt_size-1]),2);

    if (pth->compute_damping_scale == _TRUE_) {

      /* r_d scales like (1+z)**-3/2 */
      pvecthermo[pth->index_th_r_d] = pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_r_d]*pow((1+z)/(1.+pth->z_table[pth->tt_size-1]),-1.5);

    }

    /* Calculate d2kappa/dtau2 = dz/dtau d/dz[dkappa/dtau] given that [dkappa/dtau] proportional to (1+z)^2 and dz/dtau = -H */
    pvecthermo[pth->index_th_ddkappa] = -pvecback[pba->index_bg_H] * 2. / (1.+z) * pvecthermo[pth->index_th_dkappa];

    /* Calculate d3kappa/dtau3 given that [dkappa/dtau] proportional to (1+z)^2 */
    pvecthermo[pth->index_th_dddkappa] = (pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]/ (1.+z) - pvecback[pba->index_bg_H_prime]) * 2. / (1.+z) * pvecthermo[pth->index_th_dkappa];



    /* \f$ exp^{-\kappa}, g, g', g'' \f$ can be set to zero: they are
       used only for computing the source functions in the
       perturbation module; but source functions only need to be
       sampled below z_initial (the condition that
       z_start_sources<z_initial is checked in the perturbation
       module) */
    pvecthermo[pth->index_th_exp_m_kappa] = 0.;
    pvecthermo[pth->index_th_g]=0.;
    pvecthermo[pth->index_th_dg]=0.;
    pvecthermo[pth->index_th_ddg]=0.;


    /* Calculate Tb */
    pvecthermo[pth->index_th_Tb] = pba->T_cmb*(1.+z);

    /* Calculate cb2 (cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)) */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[pth->index_th_cb2] = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + x0 * (1.-pth->YHe)) * pba->T_cmb * (1.+z) * 4. / 3.;

    /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
    if (pth->compute_cb2_derivatives == _TRUE_) {

      /* since cb2 proportional to (1+z) or 1/a, its derivative wrt conformal time is given by dcb2 = - a H cb2 */
      pvecthermo[pth->index_th_dcb2] = - pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * pvecthermo[pth->index_th_cb2];

      /* then its second derivative is given by ddcb2 = - a H' cb2 */
      pvecthermo[pth->index_th_ddcb2] = - pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a] * pvecthermo[pth->index_th_cb2];
    }

    /* in this regime, variation rate = dkappa/dtau */
    pvecthermo[pth->index_th_rate] = pvecthermo[pth->index_th_dkappa];

  }

  /** - interpolate in table with array_interpolate_spline() (normal
      mode) or array_interpolate_spline_growing_closeby() (closeby
      mode) */

 else {

    /* some very specific cases require linear interpolation because of a break in the derivative of the functions */

    if ((((pth->reio_parametrization == reio_half_tanh) || (pth->reio_stars_and_dark_matter == _TRUE_))&& (z < 2*pth->z_reio))
        || ((pth->reio_parametrization == reio_inter) && (z < 50.))) {

      class_call(array_interpolate_linear(
                                          pth->z_table,
                                          pth->tt_size,
                                          pth->thermodynamics_table,
                                          pth->th_size,
                                          z,
                                          last_index,
                                          pvecthermo,
                                          pth->th_size,
                                          pth->error_message),
                 pth->error_message,
                 pth->error_message);
    }

    /* in the "normal" case, use spline interpolation */
    else {

      if (inter_mode == pth->inter_normal) {

        class_call(array_interpolate_spline(
                                            pth->z_table,
                                            pth->tt_size,
                                            pth->thermodynamics_table,
                                            pth->d2thermodynamics_dz2_table,
                                            pth->th_size,
                                            z,
                                            last_index,
                                            pvecthermo,
                                            pth->th_size,
                                            pth->error_message),
                   pth->error_message,
                   pth->error_message);
      }

      if (inter_mode == pth->inter_closeby) {

        class_call(array_interpolate_spline_growing_closeby(
                                                            pth->z_table,
                                                            pth->tt_size,
                                                            pth->thermodynamics_table,
                                                            pth->d2thermodynamics_dz2_table,
                                                            pth->th_size,
                                                            z,
                                                            last_index,
                                                            pvecthermo,
                                                            pth->th_size,
                                                            pth->error_message),
                   pth->error_message,
                   pth->error_message);

      }
    }
  }
  return _SUCCESS_;
}

/**
 * Initialize the thermo structure, and in particular the
 * thermodynamics interpolation table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input/Output: pointer to initialized thermo structure
 * @return the error status
 */
int thermodynamics_init(
                        struct precision * ppr,
                        struct background * pba,
                        struct thermo * pth
                        ) {

  /** Summary: */

  /** - define local variables */

  /* index running over time*/
  int index_tau;
  /* temporary variables related to visibility function */
  double g;
  /* vector of background values for calling background_at_tau() */
  double * pvecback;
  /* index for calling background_at_tau() */
  int last_index_back;
  /* temporary table of values of tau associated with z values in pth->z_table */
  double * tau_table;
  /* same ordered in growing time rather than growing redshift */
  double * tau_table_growing;
  /* conformal time of reionization */
  double tau_reio;
  /* R = (3./4.)*(rho_b/rho_g) */
  double R;
  /* structures for storing temporarily information on recombination
     and reionization */
  struct recombination reco;
  struct reionization reio;
  struct recombination * preco;
  struct reionization * preio;

  double tau;
  double g_max;
  int index_tau_max;

  int index_z; // GFA
  int status;
  char string1[10*_ARGUMENT_LENGTH_MAX_];
  char string2[10*_ARGUMENT_LENGTH_MAX_];
  char command_with_arguments[20*_ARGUMENT_LENGTH_MAX_];
  FILE * file_boost = fopen(ppr->boost_file, "w");



  /** - initialize pointers, allocate background vector */

  preco=&reco;
  preio=&reio;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  if (pth->thermodynamics_verbose > 0)
    printf("Computing thermodynamics");

  /** - compute and check primordial Helium fraction  */

  /* Y_He */
  if (pth->YHe == _BBN_) {
    class_call(thermodynamics_helium_from_bbn(ppr,pba,pth),
               pth->error_message,
               pth->error_message);
    if (pth->thermodynamics_verbose > 0)
      printf(" with Y_He=%.4f\n",pth->YHe);
  }
  else {
    if (pth->thermodynamics_verbose > 0)
      printf("\n");
  }

  class_test((pth->YHe < _YHE_SMALL_)||(pth->YHe > _YHE_BIG_),
             pth->error_message,
             "Y_He=%g out of bounds (%g<Y_He<%g)",pth->YHe,_YHE_SMALL_,_YHE_BIG_);

  // GFA: compute boost factor here
    if (pth->has_UCMH_spike == _TRUE_) {
        class_alloc(pth->z_table_for_boost,ppr->Number_z*sizeof(double),pth->error_message);
        class_alloc(pth->boost_table,ppr->Number_z*sizeof(double),pth->error_message);
        class_call(compute_boost_NFW_UCMH(ppr,pba,pth),pth->error_message,pth->error_message);
        /* write the results for the boost factor in a datafile, they will be read later inside DarkAges */
        for (index_z=0; index_z < ppr->Number_z; index_z++) {
         fprintf(file_boost, "%e\t %e \n",pth->z_table_for_boost[index_z], pth->boost_table[index_z]);
        }
        fclose(file_boost);
    //    strcat(ppr->command_fz," --z_boost ");
    //    for (index_z=0; index_z <= ppr->Number_z; index_z++) {
    //      sprintf(string1,"%g",pth->z_table_for_boost[index_z]);
    //      strcat(ppr->command_fz,string1);
    //      strcat(ppr->command_fz," ");
    //    }

    //    strcat(ppr->command_fz," --boost ");
    //    for (index_z=0; index_z <= ppr->Number_z; index_z++) {
    //      sprintf(string1,"%g",pth->boost_table[index_z]);
    //      strcat(ppr->command_fz,string1);
    //      strcat(ppr->command_fz," ");
    //    }

    }


  /** Initialize annihilation coefficient (First check if exotic energy injection is demanded) */

  if(pth->annihilation >0 || pth->decay_fraction > 0 || pth->PBH_accreting_mass > 0 || pth->PBH_evaporating_mass > 0 || pth->has_extended_PBH_MassFunc == _TRUE_){
    if(pth->energy_repart_coefficient==GSVI || pth->energy_repart_coefficient==no_factorization || pth->energy_repart_coefficient ==chi_from_file){

      if (pth->has_extended_PBH_MassFunc == _TRUE_) { // compute functions f(z) per channel and now also per each mass
        class_call(thermodynamics_annihilation_coefficients_init_PBH_MF(ppr,pba,pth),
                   pth->error_message,
                   pth->error_message);
      } else {
        class_call(thermodynamics_annihilation_coefficients_init(ppr,pba,pth),
                   pth->error_message,
                   pth->error_message);
      }
    }

    if(pth->has_on_the_spot==_FALSE_ && pth->energy_repart_coefficient!=no_factorization){
      class_call(thermodynamics_annihilation_f_eff_init(ppr,pba,pth,preco),
                 pth->error_message,
                 pth->error_message);
      // If we used DarkAges in the f_eff mode in the function above, there is now no conceptual
      // difference to the scenario when we just had obtained this table from a file.
      // To avoid modifying a lot of if()-clasues in the rest of the code, we adjust
      // pth->energy_deposition_function and preco->energy_deposition_function
      // as if this table was read from file.
      if (pth->energy_deposition_function == DarkAges) {
        pth->energy_deposition_function = function_from_file;
      }
    }
  }

  /** - check energy injection parameters */

  class_test((pth->annihilation<0),
	     pth->error_message,
	     "annihilation parameter cannot be negative");

// GFA: I commented out this error, to deal with p-wave DM annihilations in early mini-halos
//  class_test((pth->annihilation>1.e-4),
//	     pth->error_message,
//	     "annihilation parameter suspiciously large (%e, while typical bounds are in the range of 1e-7 to 1e-6)",
//	     pth->annihilation);

  class_test(pth->PBH_evaporating_mass > 0 && pth->PBH_evaporating_mass < 1e15 && pth->PBH_fraction > 1e-4,pth->error_message,
	     "The value of 'pth->PBH_fraction' that you enter is suspicious given the mass you chose. You are several orders of magnitude above the limit. The code doesn't handle well too high energy injection. Please choose  'pth->PBH_fraction < 1e-4'. ")


  class_test((pth->annihilation>0)&&(pba->has_cdm==_FALSE_),
        pth->error_message,
        "CDM annihilation effects require the presence of CDM!");

  // class_test((pth->annihilation_f_halo>0) && (pth->recombination==recfast),
  //            pth->error_message,
  //            "Switching on DM annihilation in halos requires using HyRec instead of RECFAST. Otherwise some values go beyond their range of validity in the RECFAST fits, and the thermodynamics module fails. Two   solutions: add 'recombination = HyRec' to your input, or set 'annihilation_f_halo = 0.' (default).");

  class_test((pth->annihilation_f_halo<0),
	     pth->error_message,
	     "Parameter for DM annihilation in halos cannot be negative");

  class_test((pth->annihilation_z_halo<0),
	     pth->error_message,
	     "Parameter for DM annihilation in halos cannot be negative");

  class_test((pth->A_spike <0), // GFA
     	 pth->error_message,
     	 "Amplitude of spike in UCMHs cannot be negative");

  class_test((pth->k_spike <0),
     	 pth->error_message,
     	 "Location of spike in UCMHs cannot be negative");

  if (pth->thermodynamics_verbose > 0)
    if ((pth->annihilation >0) && (pth->reio_parametrization == reio_none) && (ppr->recfast_Heswitch >= 3) && (pth->recombination==recfast))
      printf("Warning: if you have DM annihilation and you use recfast with option recfast_Heswitch >= 3, then the expression for CfHe_t and dy[1] becomes undefined at late times, producing nan's. This is however   masked by reionization if you are not in reio_none mode.");

  class_test((pth->decay_fraction<0),
	     pth->error_message,
	     "decay parameter cannot be negative");

  class_test((pth->decay_fraction>0)&&(pba->has_cdm==_FALSE_),
	     pth->error_message,
	     "CDM decay effects require the presence of CDM!");

  /* tests in order to prevent segmentation fault in the following */
  class_test(_not4_ == 0.,
             pth->error_message,
             "stop to avoid division by zero");
  class_test(pth->YHe == 1.,
             pth->error_message,
             "stop to avoid division by zero");
  class_test(pth->alpha_asymmetric_planck_16<1.5 || pth->alpha_asymmetric_planck_16>50 ,pth->error_message,
	     "alpha_asymmetric_planck_16 out of range [1.5,50]: rejected to avoid memory leakage.");
  /** - assign values to all indices in the structures with thermodynamics_indices()*/

  class_call(thermodynamics_indices(pth,preco,preio),
             pth->error_message,
             pth->error_message);

  /** - solve recombination and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$ with thermodynamics_recombination() */

  class_call(thermodynamics_recombination(ppr,pba,pth,preco,pvecback),
             pth->error_message,
             pth->error_message);


  /** - if there is reionization, solve reionization and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$ with thermodynamics_reionization()*/

  if ((pth->reio_parametrization != reio_none)) {
    class_call(thermodynamics_reionization(ppr,pba,pth,preco,preio,pvecback),
               pth->error_message,
               pth->error_message);
  }
  else {
    preio->rt_size=0;
    preio->index_reco_when_reio_start=-1;
  }

  /** - merge tables in recombination and reionization structures into
      a single table in thermo structure */

  class_call(thermodynamics_merge_reco_and_reio(ppr,pth,preco,preio),
             pth->error_message,
             pth->error_message);

  /** - compute table of corresponding conformal times */

  class_alloc(tau_table,pth->tt_size*sizeof(double),pth->error_message);

  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {
    class_call(background_tau_of_z(pba,
                                   pth->z_table[index_tau],
                                   tau_table+index_tau),
               pba->error_message,
               pth->error_message);
  }

  /** - store initial value of conformal time in the structure */

  pth->tau_ini = tau_table[pth->tt_size-1];

  /** - fill missing columns (quantities not computed previously but related) */

  /** - --> baryon drag interaction rate time minus one, -[1/R * kappa'], with R = 3 rho_b / 4 rho_gamma, stored temporarily in column ddkappa */

  last_index_back = pba->bg_size-1;

  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

    class_call(background_at_tau(pba,
                                 tau_table[index_tau],
                                 pba->normal_info,
                                 pba->inter_closeby,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

    R = 3./4.*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] =
      -1./R*pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa];

  }

  /** - --> second derivative of this rate, -[1/R * kappa']'', stored temporarily in column dddkappa */
  class_call(array_spline_table_line_to_line(tau_table,
                                             pth->tt_size,
                                             pth->thermodynamics_table,
                                             pth->th_size,
                                             pth->index_th_ddkappa,
                                             pth->index_th_dddkappa,
                                             _SPLINE_EST_DERIV_,
                                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> compute tau_d = [int_{tau_today}^{tau} dtau -dkappa_d/dtau] */
  class_call(array_integrate_spline_table_line_to_line(tau_table,
                                                       pth->tt_size,
                                                       pth->thermodynamics_table,
                                                       pth->th_size,
                                                       pth->index_th_ddkappa,
                                                       pth->index_th_dddkappa,
                                                       pth->index_th_tau_d,
                                                       pth->error_message),
             pth->error_message,
             pth->error_message);

  /* the temporary quantities stored in columns ddkappa and dddkappa
     will not be used anymore, they will be overwritten */

  /** - --> compute r_d = 2pi/k_d = 2pi * [int_{tau_ini}^{tau} dtau (1/kappa') (R^2+4/5(1+R))/(1+R^2)/6 ]^1/2 (see e.g. Wayne Hu's thesis eq. (5.59) */

  if (pth->compute_damping_scale == _TRUE_) {

    class_alloc(tau_table_growing,pth->tt_size*sizeof(double),pth->error_message);

    /* compute integrand 1/kappa' (R^2+4/5(1+R))/(1+R^2)/6 and store temporarily in column "ddkappa" */
    for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

      tau_table_growing[index_tau]=tau_table[pth->tt_size-1-index_tau];

      class_call(background_at_tau(pba,
                                 tau_table_growing[index_tau],
                                 pba->normal_info,
                                 pba->inter_closeby,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

      R = 3./4.*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];

      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] =
                 1./pth->thermodynamics_table[(pth->tt_size-1-index_tau)*pth->th_size+pth->index_th_dkappa]
                 *(R*R+4./5.*(1.+R))/(1.+R*R)/6.;

    }

    /* compute second derivative of integrand 1/kappa' and store temporarily in column "dddkappa" */
    class_call(array_spline_table_line_to_line(tau_table_growing,
                                               pth->tt_size,
                                               pth->thermodynamics_table,
                                               pth->th_size,
                                               pth->index_th_ddkappa,
                                               pth->index_th_dddkappa,
                                               _SPLINE_EST_DERIV_,
                                               pth->error_message),
               pth->error_message,
               pth->error_message);

    /* compute integrated quantity r_d^2 and store temporarily in column "g" */
     class_call(array_integrate_spline_table_line_to_line(tau_table_growing,
                                                          pth->tt_size,
                                                          pth->thermodynamics_table,
                                                          pth->th_size,
                                                          pth->index_th_ddkappa,
                                                          pth->index_th_dddkappa,
                                                          pth->index_th_g,
                                                          pth->error_message),
                pth->error_message,
                pth->error_message);

     free(tau_table_growing);

     /* an analytic calculation shows that in the early
        radiation-dominated and ionized universe, when kappa' is
        proportional to (1+z)^2 and tau is proportional to the scale
        factor, the integral is equal to eta/(3 kappa')*2./15. So
        [tau_ini/3/kappa'_ini*2./15.] should be added to the integral in
        order to account for the integration between 0 and tau_ini */

     /* compute r_d */
     for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_r_d] =
         2.*_PI_*sqrt(tau_table[pth->tt_size-1]/3.
                      /pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_dkappa]*2./15.
                      +pth->thermodynamics_table[(pth->tt_size-1-index_tau)*pth->th_size+pth->index_th_g]);

     }

  }

  /** - --> second derivative with respect to tau of dkappa (in view of spline interpolation) */
  class_call(array_spline_table_line_to_line(tau_table,
                                             pth->tt_size,
                                             pth->thermodynamics_table,
                                             pth->th_size,
                                             pth->index_th_dkappa,
                                             pth->index_th_dddkappa,
                                             _SPLINE_EST_DERIV_,
                                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> first derivative with respect to tau of dkappa (using spline interpolation) */
  class_call(array_derive_spline_table_line_to_line(tau_table,
                                                    pth->tt_size,
                                                    pth->thermodynamics_table,
                                                    pth->th_size,
                                                    pth->index_th_dkappa,
                                                    pth->index_th_dddkappa,
                                                    pth->index_th_ddkappa,
                                                    pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> compute -kappa = [int_{tau_today}^{tau} dtau dkappa/dtau], store temporarily in column "g" */
  class_call(array_integrate_spline_table_line_to_line(tau_table,
                                                       pth->tt_size,
                                                       pth->thermodynamics_table,
                                                       pth->th_size,
                                                       pth->index_th_dkappa,
                                                       pth->index_th_dddkappa,
                                                       pth->index_th_g,
                                                       pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives == _TRUE_) {

    /** - ---> second derivative with respect to tau of cb2 */
    class_call(array_spline_table_line_to_line(tau_table,
                                               pth->tt_size,
                                               pth->thermodynamics_table,
                                               pth->th_size,
                                               pth->index_th_cb2,
                                               pth->index_th_ddcb2,
                                               _SPLINE_EST_DERIV_,
                                               pth->error_message),
               pth->error_message,
               pth->error_message);


    /** - ---> first derivative with respect to tau of cb2 (using spline interpolation) */
    class_call(array_derive_spline_table_line_to_line(tau_table,
                                                      pth->tt_size,
                                                      pth->thermodynamics_table,
                                                      pth->th_size,
                                                      pth->index_th_cb2,
                                                      pth->index_th_ddcb2,
                                                      pth->index_th_dcb2,
                                                      pth->error_message),
               pth->error_message,
               pth->error_message);
  }


  free(tau_table);

  /** - --> compute visibility: \f$ g= (d \kappa/d \tau) e^{- \kappa} \f$ */

  /* loop on z (decreasing z, increasing time) */
  for (index_tau=pth->tt_size-1; index_tau>=0; index_tau--) {




    /** - ---> compute exp(-kappa) */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_exp_m_kappa] =
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);


    /** - ---> compute g */
    g =
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);

    /** - ---> compute g' (the plus sign of the second term is correct, see def of -kappa in thermodynamics module!) */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dg] =
      (pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] +
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa]) *
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);

    /** - ---> compute g''  */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddg] =
      (pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dddkappa] +
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] * 3. +
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa]) *
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);



    /** - ---> store g */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g] = g;

    /** - ---> compute variation rate */
    class_test(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] == 0.,
               pth->error_message,
               "variation rate diverges");

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rate] =
      sqrt(pow(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa],2)
           +pow(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa]/
                pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa],2)
           +fabs(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dddkappa]/
                 pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa]));

  }

  /** - smooth the rate (details of smoothing unimportant: only the
      order of magnitude of the rate matters) */
  class_call(array_smooth(pth->thermodynamics_table,
                          pth->th_size,
                          pth->tt_size,
                          pth->index_th_rate,
                          ppr->thermo_rate_smoothing_radius,
                          pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - fill tables of second derivatives with respect to z (in view of spline interpolation) */

  class_call(array_spline_table_lines(pth->z_table,
                                      pth->tt_size,
                                      pth->thermodynamics_table,
                                      pth->th_size,
                                      pth->d2thermodynamics_dz2_table,
                                      _SPLINE_EST_DERIV_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - find maximum of g */

  index_tau=pth->tt_size-1;
  while (pth->z_table[index_tau]>_Z_REC_MAX_) {
    index_tau--;
  }

  class_test(pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g] >
             pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g],
             pth->error_message,
             "found a recombination redshift greater or equal to the maximum value imposed in thermodynamics.h, z_rec_max=%g",_Z_REC_MAX_);

  while (pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g] <
         pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]) {
    index_tau--;
  }

  g_max = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g];
  index_tau_max = index_tau;

  /* approximation for maximum of g, using cubic interpolation, assuming equally spaced z's */
  pth->z_rec=pth->z_table[index_tau+1]+0.5*(pth->z_table[index_tau+1]-pth->z_table[index_tau])*(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]-pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g])/(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]-2.*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g]+pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g]);
  // fprintf(stdout, "z_rec %e %e %e %e %e %e\n",pth->z_table[index_tau+1],pth->z_table[index_tau],pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g],2.*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g],pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]);
  class_test((pth->z_rec+ppr->smallest_allowed_variation >= _Z_REC_MAX_),
             pth->error_message,
             "found a recombination redshift greater or equal to the maximum value imposed in thermodynamics.h, z_rec_max=%g",_Z_REC_MAX_);

  class_test((pth->z_rec-ppr->smallest_allowed_variation <= _Z_REC_MIN_),
             pth->error_message,
             "found a recombination redshift smaller or equal to the maximum value imposed in thermodynamics.h, z_rec_min=%g",_Z_REC_MIN_);

  /** - find conformal recombination time using background_tau_of_z() **/

  class_call(background_tau_of_z(pba,pth->z_rec,&(pth->tau_rec)),
             pba->error_message,
             pth->error_message);

  class_call(background_at_tau(pba,pth->tau_rec, pba->long_info, pba->inter_normal, &last_index_back, pvecback),
             pba->error_message,
             pth->error_message);

  pth->rs_rec=pvecback[pba->index_bg_rs];
  pth->ds_rec=pth->rs_rec*pba->a_today/(1.+pth->z_rec);
  pth->da_rec=pvecback[pba->index_bg_ang_distance];
  pth->ra_rec=pth->da_rec*(1.+pth->z_rec)/pba->a_today;
  pth->angular_rescaling=pth->ra_rec/(pba->conformal_age-pth->tau_rec);

  /** - find damping scale at recombination (using linear interpolation) */

  if (pth->compute_damping_scale == _TRUE_) {

    pth->rd_rec = (pth->z_table[index_tau+1]-pth->z_rec)/(pth->z_table[index_tau+1]-pth->z_table[index_tau])*pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_r_d]
      +(pth->z_rec-pth->z_table[index_tau])/(pth->z_table[index_tau+1]-pth->z_table[index_tau])*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_r_d];

  }

  /** - find time (always after recombination) at which tau_c/tau
      falls below some threshold, defining tau_free_streaming */

  class_call(background_tau_of_z(pba,pth->z_table[index_tau],&tau),
             pba->error_message,
             pth->error_message);

  while (1./pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_dkappa]/tau
         < ppr->radiation_streaming_trigger_tau_c_over_tau) {

    index_tau--;

    class_call(background_tau_of_z(pba,pth->z_table[index_tau],&tau),
               pba->error_message,
               pth->error_message);

  }

  pth->tau_free_streaming = tau;

  /** - find baryon drag time (when tau_d crosses one, using linear
      interpolation) and sound horizon at that time */

  index_tau=0;
  while ((pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_d] < 1.) && (index_tau < pth->tt_size))
    index_tau++;

  pth->z_d = pth->z_table[index_tau-1]+
    (1.-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_d])
    /(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_d]-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_d])
    *(pth->z_table[index_tau]-pth->z_table[index_tau-1]);

  class_call(background_tau_of_z(pba,pth->z_d,&(pth->tau_d)),
             pba->error_message,
             pth->error_message);

  class_call(background_at_tau(pba,pth->tau_d, pba->long_info, pba->inter_normal, &last_index_back, pvecback),
             pba->error_message,
             pth->error_message);

  pth->rs_d=pvecback[pba->index_bg_rs];
  pth->ds_d=pth->rs_d*pba->a_today/(1.+pth->z_d);

  /** - find time above which visibility falls below a given fraction of its maximum */

  index_tau=index_tau_max;
  while ((pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g] >
          g_max * ppr->neglect_CMB_sources_below_visibility)
         && (index_tau > 0))
    index_tau--;

  class_call(background_tau_of_z(pba,pth->z_table[index_tau],&(pth->tau_cut)),
             pba->error_message,
             pth->error_message);

  /** - if verbose flag set to next-to-minimum value, print the main results */

  if (pth->thermodynamics_verbose > 0) {
    printf(" -> recombination at z = %f\n",pth->z_rec);
    printf("    corresponding to conformal time = %f Mpc\n",pth->tau_rec);
    printf("    with comoving sound horizon = %f Mpc\n",pth->rs_rec);
    printf("    angular diameter distance = %f Mpc\n",pth->da_rec);
    printf("    and sound horizon angle 100*theta_s = %f\n",100.*pth->rs_rec/pth->ra_rec);
    if (pth->compute_damping_scale == _TRUE_) {
      printf("    and with comoving photon damping scale = %f Mpc\n",pth->rd_rec);
      printf("    or comoving damping wavenumber k_d = %f 1/Mpc\n",2.*_PI_/pth->rd_rec);
    }
    printf(" -> baryon drag stops at z = %f\n",pth->z_d);
    printf("    corresponding to conformal time = %f Mpc\n",pth->tau_d);
    printf("    with comoving sound horizon rs = %f Mpc\n",pth->rs_d);
    if((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
      if (pth->reio_z_or_tau==reio_tau)
        printf(" -> reionization  at z = %f\n",pth->z_reio);
      if (pth->reio_z_or_tau==reio_z)
        printf(" -> reionization with optical depth = %f\n",pth->tau_reio);
      class_call(background_tau_of_z(pba,pth->z_reio,&tau_reio),
                 pba->error_message,
                 pth->error_message);
      printf("    corresponding to conformal time = %f Mpc\n",tau_reio);
      // printf("duration of reionization = %e, with z_beg = %e, z_mid = %e, z_end = %e\n",pth->duration_of_reionization,pth->z_10_percent,pth->z_50_percent, pth->z_99_percent);
    }
    if((pth->reio_parametrization == reio_douspis_et_al) || (pth->reio_parametrization == reio_asymmetric_planck_16)){
      printf(" -> reionization with optical depth = %f\n",pth->tau_reio);
      class_call(background_tau_of_z(pba,pth->z_reio,&tau_reio),
               pba->error_message,
               pth->error_message);
      printf("    corresponding to conformal time = %f Mpc\n",tau_reio);
      // printf("duration of reionization = %e, with z_beg = %e, z_mid = %e, z_end = %e\n",pth->duration_of_reionization,pth->z_10_percent,pth->z_50_percent, pth->z_99_percent);
    }
    if (pth->reio_parametrization == reio_bins_tanh) {
      printf(" -> binned reionization gives optical depth = %f\n",pth->tau_reio);
    }
    if (pth->reio_parametrization == reio_many_tanh) {
      printf(" -> many-step reionization gives optical depth = %f\n",pth->tau_reio);
    }
    if (pth->reio_parametrization == reio_inter) {
      printf(" -> interpolated reionization history gives optical depth = %f\n",pth->tau_reio);
    }
    if (pth->thermodynamics_verbose > 1) {
      printf(" -> free-streaming approximation can be turned on as soon as tau=%g Mpc\n",
             pth->tau_free_streaming);
    }
  }

  // GFA
  if (pth->has_UCMH_spike == _TRUE_) {
    status = remove(ppr->boost_file); //remove file after all calculations are done, just to avoid generating thousands of file when running an MCMC on MontePython
    if (status == 0) {
      if (pth->thermodynamics_verbose > 1) {
        printf("Boost file deleted successfully \n");
      }
    } else {
      printf("Unable to delete the Boost file \n");
//      exit(EXIT_FAILURE);
     }
  }

  free(pvecback);

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by thermodynamics_init().
 *
 *
 * @param pth Input/Output: pointer to thermo structure (to be freed)
 * @return the error status
 */

int thermodynamics_free(
                        struct thermo * pth
                        ) {

  free(pth->z_table);
  free(pth->thermodynamics_table);
  free(pth->d2thermodynamics_dz2_table);
  if (pth->has_UCMH_spike == _TRUE_) {
    free(pth->boost_table);
    free(pth->z_table_for_boost);
  }


  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of thermodynamical quantities,
 * as well as in vector containing reionization parameters.
 *
 *
 * @param pth   Input/Output: pointer to thermo structure
 * @param preco Input/Output: pointer to recombination structure
 * @param preio Input/Output: pointer to reionization structure
 * @return the error status
 */

int thermodynamics_indices(
                           struct thermo * pth,
                           struct recombination * preco,
                           struct reionization * preio
                           ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of thermodynamics quantities */
  int index;

  /** - initialization of all indices and flags in thermo structure */
  index = 0;

  pth->index_th_xe = index;
  index++;
  pth->index_th_dkappa = index;
  index++;
  pth->index_th_tau_d = index;
  index++;
  pth->index_th_ddkappa = index;
  index++;
  pth->index_th_dddkappa = index;
  index++;
  pth->index_th_exp_m_kappa = index;
  index++;



  pth->index_th_g = index;
  index++;
  pth->index_th_dg = index;
  index++;
  pth->index_th_ddg = index;
  index++;
  pth->index_th_Tb = index;
  index++;
  pth->index_th_cb2 = index;
  index++;

  /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives == _TRUE_) {
    pth->index_th_dcb2 = index;
    index++;
    pth->index_th_ddcb2 = index;
    index++;
  }

  pth->index_th_rate = index;
  index++;

  if (pth->compute_damping_scale == _TRUE_) {
    pth->index_th_r_d = index;
    index++;
  }

  /* end of indices */
  pth->th_size = index;

  /** - initialization of all indices and flags in recombination structure */
  index = 0;

  preco->index_re_z = index;
  index++;
  preco->index_re_xe = index;
  index++;
  preco->index_re_dkappadtau = index;
  index++;
  preco->index_re_Tb = index;
  index++;
  preco->index_re_cb2 = index;
  index++;

  /* end of indices */
  preco->re_size = index;

  /** - initialization of all indices and flags in reionization structure */
  index = 0;

  preio->index_re_z = index;
  index++;
  preio->index_re_xe = index;
  index++;
  preio->index_re_Tb = index;
  index++;
  preio->index_re_cb2 = index;
  index++;
  preio->index_re_dkappadtau = index;
  index++;
  preio->index_re_dkappadz = index;
  index++;
  preio->index_re_d3kappadz3 = index;
  index++;

  /* end of indices */
  preio->re_size = index;

  /** - same with parameters of the function \f$ X_e(z)\f$ */

  index=0;

  preio->index_reio_start = index;
  index++;

  /* case where x_e(z) taken like in CAMB (other cases can be added) */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {

    preio->index_reio_redshift = index;
    index++;
    preio->index_reio_exponent = index;
    index++;
    preio->index_reio_width = index;
    index++;
    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;

    preio->reio_num_params = index;

  }

  /* case where x_e(z) is binned */
  if (pth->reio_parametrization == reio_bins_tanh) {

    /* the code will not only copy here the "bin centers" passed in
       input. It will add an initial and final value for (z,xe). So
       this array has a dimension bigger than the bin center array */

    preio->reio_num_z=pth->binned_reio_num+2; /* add two values: beginning and end of reio */

    preio->index_reio_first_z = index;
    index+= preio->reio_num_z;
    preio->index_reio_first_xe = index;
    index+= preio->reio_num_z;
    preio->index_reio_step_sharpness = index;
    index++;

  }

  /* case where x_e(z) has many tanh jumps */
  if (pth->reio_parametrization == reio_many_tanh) {

    /* the code will not only copy here the "jump centers" passed in
       input. It will add an initial and final value for (z,xe). So
       this array has a dimension bigger than the jump center array */

    preio->reio_num_z=pth->many_tanh_num+2; /* add two values: beginning and end of reio */

    preio->index_reio_first_z = index;
    index+= preio->reio_num_z;
    preio->index_reio_first_xe = index;
    index+= preio->reio_num_z;
    preio->index_reio_step_sharpness = index;
    index++;

  }

    /* case where x_e(z) must be interpolated */
  if (pth->reio_parametrization == reio_inter) {

    preio->reio_num_z=pth->reio_inter_num;

    preio->index_reio_first_z = index;
    index+= preio->reio_num_z;
    preio->index_reio_first_xe = index;
    index+= preio->reio_num_z;

  }

  if (pth->reio_parametrization == reio_douspis_et_al){

    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_lambda_douspis_et_al = index;
    index++;
    preio->index_zp_douspis_et_al = index;
    index++;
    preio->index_Qp_douspis_et_al = index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;


    preio->reio_num_params = index;

  }

  if (pth->reio_parametrization == reio_asymmetric_planck_16){

    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_alpha_asymmetric_planck_16 = index;
    index++;
    preio->index_z_end_asymmetric_planck_16= index;
    index++;
    preio->index_z_start_asymmetric_planck_16= index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;


    preio->reio_num_params = index;

  }

  if (pth->reio_parametrization == reio_stars_sfr_source_term){

    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;


    preio->reio_num_params = index;

  }

  /* flags for calling the interpolation routine */

  pth->inter_normal=0;
  pth->inter_closeby=1;

  return _SUCCESS_;
}

/**
 * Infer the primordial helium fraction from standard BBN, as a
 * function of the baryon density and expansion rate during BBN.
 *
 * This module is simpler then the one used in arXiv:0712.2826 because
 * it neglects the impact of a possible significant chemical
 * potentials for electron neutrinos. The full code with xi_nu_e could
 * be introduced here later.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input/Output: pointer to initialized thermo structure
 * @return the error status
 */
int thermodynamics_helium_from_bbn(
                                   struct precision * ppr,
                                   struct background * pba,
                                   struct thermo * pth
                                   ) {

  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;

  int num_omegab=0;
  int num_deltaN=0;

  double * omegab=NULL;
  double * deltaN=NULL;
  double * YHe=NULL;
  double * ddYHe=NULL;
  double * YHe_at_deltaN=NULL;
  double * ddYHe_at_deltaN=NULL;

  int array_line=0;
  double DeltaNeff;
  double omega_b;
  int last_index;
  double Neff_bbn, z_bbn, tau_bbn, *pvecback;

  /**Summary: */
  /** - Infer effective number of neutrinos at the time of BBN */
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  /** - 8.6173e-11 converts from Kelvin to MeV. We randomly choose 0.1 MeV to be the temperature of BBN */
  z_bbn = 0.1/(8.6173e-11*pba->T_cmb)-1.0;

  class_call(background_tau_of_z(pba,
                                 z_bbn,
                                 &tau_bbn),
             pba->error_message,
             pth->error_message);

  class_call(background_at_tau(pba,
                               tau_bbn,
                               pba->long_info,
                               pba->inter_normal,
                               &last_index,
                               pvecback),
             pba->error_message,
             pth->error_message);

  Neff_bbn = (pvecback[pba->index_bg_Omega_r]
	      *pvecback[pba->index_bg_rho_crit]
	      -pvecback[pba->index_bg_rho_g])
    /(7./8.*pow(4./11.,4./3.)*pvecback[pba->index_bg_rho_g]);

  free(pvecback);

  //  printf("Neff early = %g, Neff at bbn: %g\n",pba->Neff,Neff_bbn);




  /** - compute Delta N_eff as defined in bbn file, i.e. \f$ \Delta N_{eff}=0\f$ means \f$ N_{eff}=3.046\f$ */
  DeltaNeff = Neff_bbn - 3.046;

  /** - Since class v2.7 we make use of a fitting formula from PArthENoPE, Iocco et al. (2009), updated with the latest observational data
  on nuclear rates and neutron lifetime. It corresponds to the fit used by the Planck team, see e.g. Planck 15 cosmological papers. */

  double a[9]={0.2311,0.9502,-11.27,0.01356,0.008581,-0.1810,-0.0009795,-0.001370,0.01746}; //Coefficient of the fit.
  double b = 0.728; //Exponent of the lifetime in the fit.
  double tmp_YHe = 0.;
  int ii = 0, m = 0, n = 0;

  omega_b=pba->Omega0_b*pba->h*pba->h;

  for(ii=0;ii<9;ii++){
    if(ii == 0){
      m=0;
      n=0;
    }
    else if(ii == 3){
      m=1;
      n=0;
    }
    else if(ii==6){
      m=2;
      n=0;
    }
    tmp_YHe += pow(_NEUTRON_LIFETIME_/880.3,b)*a[ii]*pow(omega_b,n)*pow(DeltaNeff,m);
    n++;
  }

  pth->YHe = tmp_YHe;

  return _SUCCESS_;

}

/*****MODIF Vivian Poulin : Add new functions for energy repartition from DM annihilations or decays

Modification: Patrick Stöcker (20.02.17): Adding call to external script to calculate the annihilation coefficients on the fly.

*****/
int thermodynamics_annihilation_coefficients_init(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct thermo * pth
                                                  ) {

  FILE * fA = NULL;
  char line[_LINE_LENGTH_MAX_];
  char * left;

  /* BEGIN: New variables related to the use of an external code to calculate the annihilation coefficients */
  //char arguments[_ARGUMENT_LENGTH_MAX_];
  char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  int status;
  /* END */

  int num_lines=0;
  int array_line=0;

  /*

      the following file is assumed to contain (apart from comments and blank lines):
     - One number (num_lines) = number of lines of the file
     - six columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where each chi represents the fraction of energy going respectively into
     heat, excitation of lyman-alpha level, Hydrogen ionisation, Helium ionisation, photons below 10.2 eV unseeable by the IGM.

  */

  /* BEGIN: Add switch (1) */
  if (pth->energy_deposition_function == function_from_file || pth->energy_repart_coefficient == GSVI || pth->energy_repart_coefficient == chi_from_file) {
    class_open(fA,ppr->energy_injec_coeff_file, "r",pth->error_message);
  } else {
    /* Write the command */
    sprintf(command_with_arguments, "%s", ppr->command_fz);
    // free(ppr->command_fz);
    if (pth->thermodynamics_verbose > 0) {
      printf(" -> running: %s\n", command_with_arguments);
    }
    /* Launch the process and retrieve the output */
    fflush(fA); // clean the buffer
    fA = popen(command_with_arguments, "r"); // popen stands for process-open, in "read" mode
    class_test(fA == NULL, pth->error_message, "The program failed to set the environment for the external command."); //process is launched here
    // all the output of DarkAges is now stored in file fA
  }

  /* END */

  /* go through each line */
  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) { //each row of the file fA is stored in "line" variable
    /* eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* check that the line is neither blank nor a comment. In
       ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %,
       etc.) */
    if (left[0] > 39) {

      /* if the line contains data, we must interprete it. If
         num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */
      if (num_lines == 0) {

        /* read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&num_lines) != 1,
                   pth->error_message,
                   "could not read value of parameters num_lines in file %s\n",ppr->energy_injec_coeff_file); //num_lines is not 0 anymore
        class_alloc(pth->annihil_coef_xe,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_heat,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_lya,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_ionH,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_ionHe,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_lowE,num_lines*sizeof(double),pth->error_message);

        class_alloc(pth->annihil_coef_dd_heat,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_dd_lya,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_dd_ionH,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_dd_ionHe,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_dd_lowE,num_lines*sizeof(double),pth->error_message);
        pth->annihil_coef_num_lines = num_lines;


        array_line=0;

      }
      else {

        /* read coefficients */
        class_test(sscanf(line,"%lg %lg %lg %lg %lg %lg",
                          &(pth->annihil_coef_xe[array_line]),
                          &(pth->annihil_coef_heat[array_line]),
                          &(pth->annihil_coef_lya[array_line]),
                          &(pth->annihil_coef_ionH[array_line]),
                          &(pth->annihil_coef_ionHe[array_line]),
                          &(pth->annihil_coef_lowE[array_line])) != 6,
                   pth->error_message,
                   "could not read value of parameters coeeficients in file %s\n",ppr->energy_injec_coeff_file);
        if(pth->print_energy_deposition_function){
          if(array_line == 0){
                fprintf(stdout,"##################################################\n### This is the standardized output to be read by CLASS.\n### For the correct usage ensure that all other\n###'print'-commands in your script are silenced.\n##################################################\n#z_dep	f_heat	f_lya	f_ionH	f_ionHe	f_lowE\n");
          }
          printf("%e %e %e %e %e %e \n",
          (pth->annihil_coef_xe[array_line]),
          (pth->annihil_coef_heat[array_line]),
          (pth->annihil_coef_lya[array_line]),
          (pth->annihil_coef_ionH[array_line]),
          (pth->annihil_coef_ionHe[array_line]),
          (pth->annihil_coef_lowE[array_line]));
        }
        array_line ++;
      }
    }
  }
  /* BEGIN: Add switch (2) */
  if (pth->energy_deposition_function == function_from_file || pth->energy_repart_coefficient == GSVI || pth->energy_repart_coefficient == chi_from_file) {
    fclose(fA);
  } else {
    status = pclose(fA);
    class_test(status != 0., pth->error_message, "The attempt to launch the external command was not successful. Maybe the output of the external command is not in the right format.");
  }
  /* END */

  /* spline in one dimension */
  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_heat,
                                      1,
                                      pth->annihil_coef_dd_heat,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_lya,
                                      1,
                                      pth->annihil_coef_dd_lya,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_ionH,
                                      1,
                                      pth->annihil_coef_dd_ionH,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);
  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                     pth->annihil_coef_num_lines,
                                     pth->annihil_coef_ionHe,
                                     1,
                                     pth->annihil_coef_dd_ionHe,
                                     _SPLINE_NATURAL_,
                                     pth->error_message),
              pth->error_message,
              pth->error_message);
  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_lowE,
                                      1,
                                      pth->annihil_coef_dd_lowE,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
               pth->error_message,
               pth->error_message);

  return _SUCCESS_;

}

// GFA, only used in the extended PBH mass function scenarios
// this function calls DarkAges several times (one per each mass), by writting a python command with a different mass each time
// results per z and per mass are stored in the matrices annihil_coef_X_at_mass, annihil_coef_dd_X_at_mass (required for interpolation),
// where X indicates the specific channel (heat, lya, etc)
int thermodynamics_annihilation_coefficients_init_PBH_MF(
                                                         struct precision * ppr,
                                                         struct background * pba,
                                                         struct thermo * pth
                                                         ) {

  FILE * fA = NULL;
  char line[_LINE_LENGTH_MAX_];
  char * left;

  char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  char string2[_ARGUMENT_LENGTH_MAX_];
  int status;
  int index_M;

  int num_lines;
  int array_line;

  class_alloc(pth->annihil_coef_xe_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);
  class_alloc(pth->annihil_coef_heat_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);
  class_alloc(pth->annihil_coef_lya_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);
  class_alloc(pth->annihil_coef_ionH_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);
  class_alloc(pth->annihil_coef_ionHe_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);
  class_alloc(pth->annihil_coef_lowE_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);

  class_alloc(pth->annihil_coef_dd_heat_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);
  class_alloc(pth->annihil_coef_dd_lya_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);
  class_alloc(pth->annihil_coef_dd_ionH_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);
  class_alloc(pth->annihil_coef_dd_ionHe_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);
  class_alloc(pth->annihil_coef_dd_lowE_at_mass, sizeof(double*)*pth->num_PBH_accreting_mass,pba->error_message);

  class_alloc(pth->annihil_coef_num_lines_at_mass, sizeof(int)*pth->num_PBH_accreting_mass,pba->error_message);


  for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) {

    num_lines= 0;
    array_line=0;
    sprintf(ppr->command_fz,"");
    strcat(ppr->command_fz, "python ");
    strcat(ppr->command_fz,__CLASSDIR__);
    strcat(ppr->command_fz,"/DarkAgesModule/bin/DarkAges --hist=accreting_PBH --mass=");
    sprintf(string2,"%g",pth->table_PBH_accreting_mass[index_M]);
    strcat(ppr->command_fz,string2);
    if(pth->PBH_accretion_recipe==spherical_accretion) {
      strcat(ppr->command_fz," --accretion_recipe=spherical_accretion");
    } else if (pth->PBH_accretion_recipe==disk_accretion) {
      strcat(ppr->command_fz," --accretion_recipe=disk_accretion");
    }
    strcat(ppr->command_fz," --Log10Emin=0 --Log10Emax=5.5 --nbins_table=20");
    sprintf(command_with_arguments, "%s", ppr->command_fz);
    if (pth->thermodynamics_verbose > 0) {
      printf(" -> running: %s\n", command_with_arguments);
    }
    fflush(fA);
    fA = popen(command_with_arguments, "r");
    class_test(fA == NULL, pth->error_message, "The program failed to set the environment for the external command.");

    while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
      left=line;
      while (left[0]==' ') {
        left++;
      }

      if (left[0] > 39) {
        if (num_lines == 0) {
          class_test(sscanf(line,"%d",&num_lines) != 1,
                     pth->error_message,
                     "could not read value of parameters num_lines in file %s\n",ppr->energy_injec_coeff_file); //num_lines is not 0 anymore
          class_alloc(pth->annihil_coef_xe_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          class_alloc(pth->annihil_coef_heat_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          class_alloc(pth->annihil_coef_lya_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          class_alloc(pth->annihil_coef_ionH_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          class_alloc(pth->annihil_coef_ionHe_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          class_alloc(pth->annihil_coef_lowE_at_mass[index_M],num_lines*sizeof(double),pth->error_message);

          class_alloc(pth->annihil_coef_dd_heat_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          class_alloc(pth->annihil_coef_dd_lya_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          class_alloc(pth->annihil_coef_dd_ionH_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          class_alloc(pth->annihil_coef_dd_ionHe_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          class_alloc(pth->annihil_coef_dd_lowE_at_mass[index_M],num_lines*sizeof(double),pth->error_message);
          pth->annihil_coef_num_lines_at_mass[index_M] = num_lines;
          array_line=0;
        }
        else {
          /* read coefficients */
          class_test(sscanf(line,"%lg %lg %lg %lg %lg %lg",
                            &(pth->annihil_coef_xe_at_mass[index_M][array_line]),
                            &(pth->annihil_coef_heat_at_mass[index_M][array_line]),
                            &(pth->annihil_coef_lya_at_mass[index_M][array_line]),
                            &(pth->annihil_coef_ionH_at_mass[index_M][array_line]),
                            &(pth->annihil_coef_ionHe_at_mass[index_M][array_line]),
                            &(pth->annihil_coef_lowE_at_mass[index_M][array_line])) != 6,
                     pth->error_message,
                     "could not read value of parameters coeeficients in file \n");
          array_line ++;
        }
      }
    }

    status = pclose(fA);
    class_test(status != 0., pth->error_message, "The attempt to launch the external command was not successful. Maybe the output of the external command is not in the right format.");

    class_call(array_spline_table_lines(pth->annihil_coef_xe_at_mass[index_M],
                                        pth->annihil_coef_num_lines_at_mass[index_M],
                                        pth->annihil_coef_heat_at_mass[index_M],
                                        1,
                                        pth->annihil_coef_dd_heat_at_mass[index_M],
                                        _SPLINE_NATURAL_,
                                        pth->error_message),
               pth->error_message,
               pth->error_message);

    class_call(array_spline_table_lines(pth->annihil_coef_xe_at_mass[index_M],
                                        pth->annihil_coef_num_lines_at_mass[index_M],
                                        pth->annihil_coef_lya_at_mass[index_M],
                                        1,
                                        pth->annihil_coef_dd_lya_at_mass[index_M],
                                        _SPLINE_NATURAL_,
                                        pth->error_message),
               pth->error_message,
               pth->error_message);

    class_call(array_spline_table_lines(pth->annihil_coef_xe_at_mass[index_M],
                                        pth->annihil_coef_num_lines_at_mass[index_M],
                                        pth->annihil_coef_ionH_at_mass[index_M],
                                        1,
                                        pth->annihil_coef_dd_ionH_at_mass[index_M],
                                        _SPLINE_NATURAL_,
                                        pth->error_message),
               pth->error_message,
               pth->error_message);
    class_call(array_spline_table_lines(pth->annihil_coef_xe_at_mass[index_M],
                                       pth->annihil_coef_num_lines_at_mass[index_M],
                                       pth->annihil_coef_ionHe_at_mass[index_M],
                                       1,
                                       pth->annihil_coef_dd_ionHe_at_mass[index_M],
                                       _SPLINE_NATURAL_,
                                       pth->error_message),
                pth->error_message,
                pth->error_message);
    class_call(array_spline_table_lines(pth->annihil_coef_xe_at_mass[index_M],
                                        pth->annihil_coef_num_lines_at_mass[index_M],
                                        pth->annihil_coef_lowE_at_mass[index_M],
                                        1,
                                        pth->annihil_coef_dd_lowE_at_mass[index_M],
                                        _SPLINE_NATURAL_,
                                        pth->error_message),
                 pth->error_message,
                 pth->error_message);

  }

  return _SUCCESS_;
}


/**
 * This function is used by the energy injection module for two different interpolations:
 * Either to directly interpolate the f(z) functions per channels or to interpolate
 * the chi(x_e) functions  per channels when the factorisation approximation is assumed.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input/Output: pointer to initialized thermo structure
 * @return the error status
 */

int thermodynamics_annihilation_coefficients_interpolate(
                                                         struct precision * ppr,
                                                         struct background * pba,
                                                         struct thermo * pth,
                                                         double xe_or_z
                                                         ) {

  int last_index;

  class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_heat,
                                      pth->annihil_coef_dd_heat,
                                      1,
                                      xe_or_z,
                                      &last_index,
                                      &(pth->chi_heat),
                                      1,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_lya,
                                      pth->annihil_coef_dd_lya,
                                      1,
                                      xe_or_z,
                                      &last_index,
                                      &(pth->chi_lya),
                                      1,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_ionH,
                                      pth->annihil_coef_dd_ionH,
                                      1,
                                      xe_or_z,
                                      &last_index,
                                      &(pth->chi_ionH),
                                      1,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);
    class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                       pth->annihil_coef_num_lines,
                                       pth->annihil_coef_ionHe,
                                       pth->annihil_coef_dd_ionHe,
                                       1,
                                       xe_or_z,
                                       &last_index,
                                       &(pth->chi_ionHe),
                                       1,
                                       pth->error_message),
              pth->error_message,
              pth->error_message);
    class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                       pth->annihil_coef_num_lines,
                                       pth->annihil_coef_lowE,
                                       pth->annihil_coef_dd_lowE,
                                       1,
                                       xe_or_z,
                                       &last_index,
                                       &(pth->chi_lowE),
                                       1,
                                       pth->error_message),
              pth->error_message,
              pth->error_message);

  return _SUCCESS_;

}

// GFA, in the extended PBH mass function scenario this interpolate the f(z) functions per channel,
// for each of the different masses considered
int thermodynamics_annihilation_coefficients_interpolate_PBH_MF(
                                                                struct precision * ppr,
                                                                struct background * pba,
                                                                struct thermo * pth,
                                                                double xe_or_z
                                                                ) {

  int last_index;
  int index_M;

  for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) {

    class_call(array_interpolate_spline(pth->annihil_coef_xe_at_mass[index_M],
                                        pth->annihil_coef_num_lines_at_mass[index_M],
                                        pth->annihil_coef_heat_at_mass[index_M],
                                        pth->annihil_coef_dd_heat_at_mass[index_M],
                                        1,
                                        xe_or_z,
                                        &last_index,
                                        &(pth->chi_heat_at_mass[index_M]),
                                        1,
                                        pth->error_message),
               pth->error_message,
               pth->error_message);

    class_call(array_interpolate_spline(pth->annihil_coef_xe_at_mass[index_M],
                                        pth->annihil_coef_num_lines_at_mass[index_M],
                                        pth->annihil_coef_lya_at_mass[index_M],
                                        pth->annihil_coef_dd_lya_at_mass[index_M],
                                        1,
                                        xe_or_z,
                                        &last_index,
                                        &(pth->chi_lya_at_mass[index_M]),
                                        1,
                                        pth->error_message),
               pth->error_message,
               pth->error_message);

    class_call(array_interpolate_spline(pth->annihil_coef_xe_at_mass[index_M],
                                        pth->annihil_coef_num_lines_at_mass[index_M],
                                        pth->annihil_coef_ionH_at_mass[index_M],
                                        pth->annihil_coef_dd_ionH_at_mass[index_M],
                                        1,
                                        xe_or_z,
                                        &last_index,
                                        &(pth->chi_ionH_at_mass[index_M]),
                                        1,
                                        pth->error_message),
               pth->error_message,
               pth->error_message);
      class_call(array_interpolate_spline(pth->annihil_coef_xe_at_mass[index_M],
                                         pth->annihil_coef_num_lines_at_mass[index_M],
                                         pth->annihil_coef_ionHe_at_mass[index_M],
                                         pth->annihil_coef_dd_ionHe_at_mass[index_M],
                                         1,
                                         xe_or_z,
                                         &last_index,
                                         &(pth->chi_ionHe_at_mass[index_M]),
                                         1,
                                         pth->error_message),
                pth->error_message,
                pth->error_message);
      class_call(array_interpolate_spline(pth->annihil_coef_xe_at_mass[index_M],
                                         pth->annihil_coef_num_lines_at_mass[index_M],
                                         pth->annihil_coef_lowE_at_mass[index_M],
                                         pth->annihil_coef_dd_lowE_at_mass[index_M],
                                         1,
                                         xe_or_z,
                                         &last_index,
                                         &(pth->chi_lowE_at_mass[index_M]),
                                         1,
                                         pth->error_message),
                pth->error_message,
                pth->error_message);

  }

  return _SUCCESS_;

}


int thermodynamics_annihilation_coefficients_free(
                                                  struct thermo * pth
                                                  ) {

  free(pth->annihil_coef_xe);
  free(pth->annihil_coef_heat);
  free(pth->annihil_coef_lya);
  free(pth->annihil_coef_ionH);
  free(pth->annihil_coef_ionHe);
  free(pth->annihil_coef_lowE);

  free(pth->annihil_coef_dd_heat);
  free(pth->annihil_coef_dd_lya);
  free(pth->annihil_coef_dd_ionH);
  free(pth->annihil_coef_dd_ionHe);
  free(pth->annihil_coef_dd_lowE);


  return _SUCCESS_;

}
// /****************MODIF Vivian Poulin 2 : Add f_eff function****************/
int thermodynamics_annihilation_f_eff_init(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct thermo * pth,
                                                  struct recombination * preco
                                                  ) {

  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;

  int num_lines=0;
  int array_line=0;

  /*

      the following file (or output of DarkAges) is assumed to contain (apart from comments and blank lines):
     - One number (num_lines) = number of lines of the file
     - One column (z , f(z)) where f(z) represents the "effective" fraction of energy deposited into the medium at redshift z, in presence of halo formation.

  */
  if (pth->energy_deposition_function==DarkAges) {
    if (pth->thermodynamics_verbose > 0) {
      printf(" -> running: %s\n", ppr->command_fz);
    }
    fflush(fA);
    fA = popen(ppr->command_fz, "r");
    class_test(fA == NULL, pth->error_message, "The program failed to set the environment for the external command.");
  } else {
    class_open(fA,ppr->energy_injec_f_eff_file, "r",pth->error_message);
  }

  /* go through each line */
  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {

    /* eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* check that the line is neither blank nor a comment. In
       ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %,
       etc.) */
    if (left[0] > 39) {

      /* if the line contains data, we must interprete it. If
         num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */
      if (num_lines == 0) {

        /* read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&num_lines) != 1,
                   pth->error_message,
                   "could not read value of parameters num_lines in file %s\n",ppr->energy_injec_f_eff_file);
        class_alloc(preco->annihil_z,num_lines*sizeof(double),pth->error_message);
        class_alloc(preco->annihil_f_eff,num_lines*sizeof(double),pth->error_message);

        class_alloc(preco->annihil_dd_f_eff,num_lines*sizeof(double),pth->error_message);

        preco->annihil_f_eff_num_lines = num_lines;


        array_line=0;

      }
      else {

        /* read coefficients */
        class_test(sscanf(line,"%lg %lg",
                          &(preco->annihil_z[array_line]),
                          &(preco->annihil_f_eff[array_line]))!= 2,
                   pth->error_message,
                   "could not read value of parameters coefficients in file %s\n",ppr->energy_injec_f_eff_file);
        array_line ++;
      }
    }
  }

  if (pth->energy_deposition_function == DarkAges) {
    pclose(fA);
  } else {
    fclose(fA);
  }

  /* spline in one dimension */
  class_call(array_spline_table_lines(preco->annihil_z,
                                      preco->annihil_f_eff_num_lines,
                                      preco->annihil_f_eff,
                                      1,
                                      preco->annihil_dd_f_eff,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  return _SUCCESS_;

}


int thermodynamics_annihilation_f_eff_interpolate(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z,
                                                  ErrorMsg error_message
                                                ) {

  int last_index;
  class_call(array_interpolate_spline(preco->annihil_z,
                                      preco->annihil_f_eff_num_lines,
                                      preco->annihil_f_eff,
                                      preco->annihil_dd_f_eff,
                                      1,
                                      z,
                                      &last_index,
                                      &(preco->f_eff),
                                      1,
                                      error_message),
             error_message,
             error_message);


  return _SUCCESS_;

}


int thermodynamics_annihilation_f_eff_free(
                                                  struct recombination * preco
                                                  ) {

  free(preco->annihil_z);
  free(preco->annihil_f_eff);
  free(preco->annihil_dd_f_eff);


  return _SUCCESS_;

}
/********** END OF MODIFICATION By Vivian Poulin **************/
/**
 * In case of non-minimal cosmology, this function determines the
 * energy rate injected in the IGM at a given redshift z (= on-the-spot
 * annihilation). This energy injection may come e.g. from dark matter
 * annihilation or decay.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param preco Input: pointer to recombination structure
 * @param z Input: redshift
 * @param energy_rate Output: energy density injection rate
 * @param error_message Output: error message
 * @return the error status
 */

/**********************************************************************************************/
/******************************Energy Injection DM annihilation**********************************/
int thermodynamics_DM_annihilation_energy_injection(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z,
                                                  double * energy_rate,
                                                  ErrorMsg error_message
                                                ){

  double rho_cdm_today;
  double Boost_factor;

  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */

  if (preco->has_UCMH_spike == _TRUE_) {
    if (z>1.e-3) {
      Boost_factor = array_interpolate_linear_simpler(preco->z_table_for_boost,ppr->Number_z,preco->boost_table,z);
    } else {
      Boost_factor = preco->boost_table[0];
    }

    *energy_rate = pow(rho_cdm_today,2)/_c_/_c_*(pow((1.+z),6)*preco->annihilation)*Boost_factor;
    /* energy density rate in J/m^3/s (remember that sigma_thermal/(preco->annihilation_m_DM*conversion) is in m^3/s/Kg) */

  } else {

    if(preco->annihilation_z_halo>0.) {
      Boost_factor = preco->annihilation_f_halo*erfc((1+z)/(1+preco->annihilation_z_halo))/pow(1+z,3);
    }
     else Boost_factor = 0;

    *energy_rate = pow(rho_cdm_today,2)/_c_/_c_*(pow((1.+z),6)*preco->annihilation)*(1+Boost_factor);
    /* energy density rate in J/m^3/s (remember that sigma_thermal/(preco->annihilation_m_DM*conversion) is in m^3/s/Kg) */

  }

}
/******************************Energy Injection DM decay**********************************/
int thermodynamics_DM_decay_energy_injection(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z,
                                                  double * energy_rate,
                                                  ErrorMsg error_message
                                                ){
  double rho_cdm_today, rho_dcdm,decay_factor;
  double tau;
  int last_index_back;
  double * pvecback;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  if(pba->Omega_ini_dcdm!=0 || pba->Omega0_dcdmdr !=0){
    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               ppr->error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               ppr->error_message);
     rho_dcdm = pvecback[pba->index_bg_rho_dcdm]*pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in J/m^3 */
     /* If uncommented, these lines allow to check approximation when computing the dcdm density with analytical results. Works very well until Omega_lambda dominates, then ~10% difference. */
    //  result_integrale = exp(-pba->Gamma_dcdm*2*((pba->Omega0_b+pba->Omega0_cdm)*pow(pba->Omega0_g+(pba->Omega0_b+pba->Omega0_cdm)/(1+z),0.5)
    //  +2*pow(pba->Omega0_g,1.5)*(1+z)-2*pba->Omega0_g*pow((1+z)*(pba->Omega0_g*(1+z)+(pba->Omega0_b+pba->Omega0_cdm)),0.5))/(3*pow((pba->Omega0_b+pba->Omega0_cdm),2)*(1+z)*pba->H0));
    //   rho_dcdm_approchee = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega_ini_dcdm)*_c_*_c_*result_integrale*pow(1+z,3);
    //   fprintf(stdout, "z = %e vrai = %e  approchee = %e relativ diff = %e\n",z,rho_dcdm, rho_dcdm_approchee,(rho_dcdm-rho_dcdm_approchee)/rho_dcdm_approchee);
  }
  else{
    if(preco->has_on_the_spot == _FALSE_)decay_factor=1; //The effect of the exponential decay is already incorporated within the f_z functions.
    else {
    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               ppr->error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               ppr->error_message);
    decay_factor = exp(-pba->Gamma_dcdm*pvecback[pba->index_bg_time]);
    }
    rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */
    rho_dcdm = rho_cdm_today*pow((1+z),3)*decay_factor; // This trick avoid mixing gravitational and electromagnetic impacts of the decay on the CMB power spectra.
  }


  *energy_rate = rho_dcdm*preco->decay_fraction*(pba->Gamma_dcdm*_c_/_Mpc_over_m_);
  // fprintf(stdout, "*energy_rate %e\n",*energy_rate );
  free(pvecback);


}
int PBH_evaporating_mass_time_evolution(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct recombination * preco,
                                  ErrorMsg error_message
                                ){

        double loop_z, loop_tau, current_mass, time_now, time_prev, dt, dz, f;
        double QCD_activation, current_pbh_temperature;
        double * pvecback_loop;
        int  i_step, last_index_back_loop;
        time_prev = 0.;
        preco->PBH_z_evaporation = 0;
        preco->PBH_table_size = ppr->recfast_Nz0;
        dz = ppr->recfast_z_initial /(preco->PBH_table_size);
        loop_z = ppr->recfast_z_initial-dz;
        current_mass = preco->PBH_evaporating_mass;
        class_alloc(pvecback_loop,pba->bg_size*sizeof(double),pba->error_message);
        class_alloc(preco->PBH_table_z,preco->PBH_table_size*sizeof(double),error_message);
        class_alloc(preco->PBH_table_mass,preco->PBH_table_size*sizeof(double),error_message);
        class_alloc(preco->PBH_table_mass_dd,preco->PBH_table_size*sizeof(double),error_message);
        class_alloc(preco->PBH_table_F,preco->PBH_table_size*sizeof(double),error_message);
        class_alloc(preco->PBH_table_F_dd,preco->PBH_table_size*sizeof(double),error_message);
        for(i_step = 0; i_step < preco->PBH_table_size; i_step++) {
          /*
          For the parametrization of F(M) we follow PRD44 (1991) 376 with the additional
          modification that we dress the "free QCD-particles" (gluons and quarks)
          with an sigmoid-activation function
          (in log10-space: Mean at 0.3 GeV and a with of 0.1*"order of magnitude")
          and the hadrons with 1 - activation to take the QCD-phase transition into account
          and to be in agreement with PRD41 (1990) 3052, where the Ansatz is taken that
          a black hole emmits those particles which appear elementary at the given energy.

          The order of the particles in the following definition of f:
          photon, neutrino, electron, muon, tau, up, down, charm, strange, top, bottom, W, Z, gluon, Higgs, neutral Pion and charged pion
          */
          current_pbh_temperature = 1.06e13 / current_mass;
          QCD_activation = 1 / (1 + exp( -(log(current_pbh_temperature)-log(0.3))/(log(10)*0.1) ));
          f =  2 * 0.060 \
            +  6 * 0.147 \
            +  4 * 0.142 * exp(-(current_mass * 5.11e-4) / (4.53 * 1.06e13)) \
            +  4 * 0.142 * exp(-(current_mass * 0.1037)  / (4.53 * 1.06e13)) \
            +  4 * 0.142 * exp(-(current_mass *  1.777)  / (4.53 * 1.06e13)) \
            + 12 * 0.142 * exp(-(current_mass * 2.2e-3)  / (4.53 * 1.06e13)) * QCD_activation \
            + 12 * 0.142 * exp(-(current_mass * 4.7e-3)  / (4.53 * 1.06e13)) * QCD_activation \
            + 12 * 0.142 * exp(-(current_mass * 1.82)    / (4.53 * 1.06e13)) * QCD_activation \
            + 12 * 0.142 * exp(-(current_mass * 9.6e-2)  / (4.53 * 1.06e13)) * QCD_activation \
            + 12 * 0.142 * exp(-(current_mass * 173.1)   / (4.53 * 1.06e13)) * QCD_activation \
            + 12 * 0.142 * exp(-(current_mass * 4.18)    / (4.53 * 1.06e13)) * QCD_activation \
            +  6 * 0.060 * exp(-(current_mass * 80.39)   / (6.04 * 1.06e13)) \
            +  3 * 0.060 * exp(-(current_mass * 91.19)   / (6.04 * 1.06e13)) \
            + 16 * 0.060 * exp(-(current_mass * 6e-1)    / (6.04 * 1.06e13)) * QCD_activation \
            +  1 * 0.267 * exp(-(current_mass * 125.06)  / (2.66 * 1.06e13)) \
            +  1 * 0.267 * exp(-(current_mass * 1.350e-1)/ (2.66 * 1.06e13)) * (1 - QCD_activation) \
            +  2 * 0.267 * exp(-(current_mass * 1.396e-1)/ (2.66 * 1.06e13)) * (1 - QCD_activation);

          class_call(background_tau_of_z(pba,
                 loop_z,
                 &loop_tau),
         pba->error_message,
         ppr->error_message);
          class_call(background_at_tau(pba,
               loop_tau,
               pba->long_info,
               pba->inter_normal,
               &last_index_back_loop,
               pvecback_loop),
         pba->error_message,
         ppr->error_message);
          time_now = pvecback_loop[pba->index_bg_time]/(_c_ / _Mpc_over_m_);
          dt = time_now - time_prev;
          time_prev = time_now;
          if (i_step > 0) {
            if (current_mass > 0.5*preco->PBH_evaporating_mass) {
              current_mass = current_mass - 5.34e-5*f*pow(current_mass/1e10,-2)*1e10 * dt;
            }
            else {
              if(preco->PBH_z_evaporation == 0)preco->PBH_z_evaporation=loop_z;
              current_mass = 0.;
              f = 0.;
            }
          }

          preco->PBH_table_z[i_step] = loop_z;
          preco->PBH_table_mass[i_step] = current_mass;
          preco->PBH_table_F[i_step] = f;
          loop_z = MAX(0,loop_z-dz);
          // printf("f %e mass %e \n",f,current_mass);
        }
        free(pvecback_loop);
        class_call(array_spline_table_lines(preco->PBH_table_z,
                    preco->PBH_table_size,
              preco->PBH_table_mass,
              1,
              preco->PBH_table_mass_dd,
              _SPLINE_NATURAL_,
              error_message),
             error_message,
             error_message);
        class_call(array_spline_table_lines(preco->PBH_table_z,
                    preco->PBH_table_size,
              preco->PBH_table_F,
              1,
              preco->PBH_table_F_dd,
              _SPLINE_NATURAL_,
              error_message),
             error_message,
             error_message);

}
/******************************Energy Injection low mass PBH (evaporation)**********************************/
int thermodynamics_evaporating_pbh_energy_injection(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z,
                                                  double * energy_rate,
                                                  ErrorMsg error_message
                                                ){

  double rho_cdm_today;
  //double tau;
  int last_index_back;
  //Parameters related to PBH

  double f, f_neutrinos, em_branching, pbh_mass;
  double dMdt;


  /* Calculate the PBH-mass evolution at first call of the function */
  if ((preco->PBH_table_is_initialized) == _FALSE_) {
    preco->PBH_table_is_initialized = _TRUE_;
    PBH_evaporating_mass_time_evolution(ppr,pba,preco,error_message);
  }
  /* End of PBH-mass loop */

  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */

  class_test(preco->PBH_table_is_initialized == _FALSE_, error_message, "The PBH table is not initialized");
  class_call(array_interpolate_spline(preco->PBH_table_z,
				      preco->PBH_table_size,
				      preco->PBH_table_mass,
				      preco->PBH_table_mass_dd,
				      1,
				      z,
				      &last_index_back,
				      &(pbh_mass),
				      1,
				      error_message),
	     error_message,
	     error_message);
  class_call(array_interpolate_spline(preco->PBH_table_z,
				      preco->PBH_table_size,
				      preco->PBH_table_F,
				      preco->PBH_table_F_dd,
				      1,
				      z,
				      &last_index_back,
				      &f,
				      1,
				      error_message),
	     error_message,
	     error_message);
  //  printf("preco: z %e pbhmass %e f %e\n",z,pbh_mass,f);

  f_neutrinos = 6*0.147;
  // em_branching = (f-f_neutrinos)/f;
  em_branching = 1.; // Currently incoporated in the computation of the f(z) functions.
  // printf("preco->PBH_z_evaporation %e\n", preco->PBH_z_evaporation);
  if(pbh_mass <= 0.0001*preco->PBH_evaporating_mass || f <= 0 || isnan(pbh_mass)==1 || isnan(f)==1 || z < preco->PBH_z_evaporation){
    pbh_mass = 0;
    dMdt = 0;
    f = 0.;
  }
  else {
    dMdt=5.34e-5*f*pow(pbh_mass/1e10,-2)*1e10;
  }
  *energy_rate = rho_cdm_today*pow((1+z),3)*preco->PBH_fraction/preco->PBH_evaporating_mass*em_branching*(dMdt);

  if(isnan(*energy_rate)==1 || *energy_rate < 0){
    *energy_rate=0.;
  }
  // if(pbh_mass>0)fprintf(stdout,"z = %lg | f = %lg | mass = %lg | energy_rate = %lg\n",z,f,pbh_mass,*energy_rate);
  // if(pbh_mass>0)fprintf(stdout,"%e %e %e %e \n",z,f,pbh_mass,*energy_rate);
}
/******************************Energy Injection high mass PBH (accretion)**********************************/
int thermodynamics_accreting_pbh_energy_injection(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z,
                                                  double * energy_rate,
                                                  ErrorMsg error_message
                                                ){

  double rho_cdm_today;
  double tau;
  int last_index_back;
  double * pvecback;
  //Parameters related to PBH
  double c_s, v_eff,v_eff_2,v_l, r_B,x_e,beta,beta_eff,beta_hat,x_cr,lambda,n_gas,M_b_dot,M_sun,M_ed_dot,epsilon,L_acc,Integrale,Normalization;
  double m_H, m_dot, m_dot_2, L_acc_2,L_ed,l,l2,M_crit;
  double rho, m_p = 938, m_e = 0.511, T_infinity = 0, rho_infinity = 0, x_e_infinity = 0, P_infinity = 0, rho_cmb = 0, t_B = 0, v_B = 0;
  double lambda_1,lambda_2,lambda_ad,lambda_iso,gamma_cooling,beta_compton_drag, T_s, T_ion, Y_s, J,tau_cooling;
  double Value_min, Value_med, Value_max, a=0, epsilon_0=0.1;
  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */

  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             error_message);


        c_s = 5.7e3*pow(preco->Tm_tmp/2730,0.5);//conversion km en m
        M_sun = 2e30; // in Kg
        n_gas = 200*1e6*pow((1+z)/1000,3); // 1e6 = conversion cm^-3 en m^-3;
        m_H= 1.67e-27; // Hydrogen mass in kg

        // x_e = 1;
        x_e = preco->xe_tmp;
        T_infinity = preco->Tm_tmp*_eV_over_Kelvin_*1e-6; //Temperature in MeV
        /** Disk accretion from Poulin et al. 1707.04206 */
        if(preco->PBH_accretion_recipe == disk_accretion){
            L_ed = 4*_PI_*_G_*preco->PBH_accreting_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);
            M_ed_dot= 10*L_ed/(_c_*_c_);
            M_crit = 0.01*M_ed_dot;
            v_B = sqrt((1+x_e)*T_infinity/m_p)*_c_;
            if(preco->PBH_relative_velocities < 0.){
              v_l = 30*MIN(1,z/1000)*1e3; // in m/s.
              if(v_B < v_l) v_eff = sqrt(v_B*v_l);
              else v_eff = v_B;
            }
            else{
              v_l = preco->PBH_relative_velocities*1e3; // converted to m/s.
              v_eff = pow(v_l*v_l+v_B*v_B,0.5);
            }

            lambda = preco->PBH_accretion_eigenvalue;
            rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in kg/m^3 */
            M_b_dot = 4*_PI_*lambda*pow(_G_*preco->PBH_accreting_mass*M_sun,2)*rho*pow(v_eff,-3.);
            if(preco->PBH_ADAF_delta == 1e-3){
              Value_min = 7.6e-5;
              Value_med = 4.5e-3;
              Value_max = 7.1e-3;
              if(M_b_dot/M_ed_dot <= Value_min){
                epsilon_0 = 0.065;
                a = 0.71;
              }
              else if(Value_min < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot  <= Value_med){
                epsilon_0 = 0.020;
                a = 0.47;
              }
              else if(Value_med < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot <= Value_max){
                epsilon_0 = 0.26;
                a = 3.67;
              }
              else{
                epsilon_0 = 0.1;
                a = 0.;
              }
              epsilon = epsilon_0 * pow(M_b_dot / M_crit,a);
            }
            else if (preco->PBH_ADAF_delta == 0.1){
              Value_min = 9.4e-5;
              Value_med = 5e-3;
              Value_max = 6.6e-3;
              if(M_b_dot/M_ed_dot <= Value_min){
                epsilon_0 = 0.12;
                a = 0.59;
              }
              else if(Value_min < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot  <= Value_med){
                epsilon_0 = 0.026;
                a = 0.27;
              }
              else if(Value_med < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot <= Value_max){
                epsilon_0 = 0.50;
                a = 4.53;
              }
              else{
                epsilon_0 = 0.1;
                a = 0.;
              }
              epsilon = epsilon_0 * pow(M_b_dot / M_crit,a);
            }
            else if (preco->PBH_ADAF_delta == 0.5){

              Value_min = 2.9e-5;
              Value_med = 3.3e-3;
              Value_max = 5.3e-3;
              if(M_b_dot/M_ed_dot <= Value_min){
                epsilon_0 = 1.58;
                a = 0.65;
              }
              else if(Value_min < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot  <= Value_med){
                epsilon_0 = 0.055;
                a = 0.076;
              }
              else if(Value_med < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot <= Value_max){
                epsilon_0 = 0.17;
                a = 1.12;
              }
              else{
                epsilon_0 = 0.1;
                a = 0.;
              }
              epsilon = epsilon_0 * pow(M_b_dot / M_crit,a);
            }

            L_acc = epsilon*M_b_dot*_c_*_c_;

          }
        /** Spherical accretion from Ali-Haimoud et al. 1612.05644 */
        else if(preco->PBH_accretion_recipe == spherical_accretion){
          rho_cmb = pvecback[pba->index_bg_rho_g]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_*_c_*_c_* 6.241509e12; /* energy density in MeV/m^3 */
          // x_e_infinity = 1; // change to 1 for the strong-feedback case
          x_e_infinity = x_e; // change to x_e for the no-feedback case
          v_B = sqrt((1+x_e_infinity)*T_infinity/m_p)*_c_; //sound speed.
          if(preco->PBH_relative_velocities < 0.){
            v_l = 30*MIN(1,z/1000)*1e3; // in m/s.
            if(v_B < v_l) v_eff = sqrt(v_B*v_l);
            else v_eff = v_B;
          }
          else{
            v_l = preco->PBH_relative_velocities*1e3; // converted to m/s.
            v_eff = pow(v_l*v_l+v_B*v_B,0.5);
          }
          r_B = _G_*preco->PBH_accreting_mass*M_sun*pow(v_eff,-2); // in m
          t_B = _G_*preco->PBH_accreting_mass*M_sun/pow(v_eff,3); // in s
          beta_compton_drag = 4./3*x_e_infinity*_sigma_*rho_cmb*t_B/(m_p)*_c_;
          gamma_cooling = 2*m_p/(m_e*(1+x_e_infinity))*beta_compton_drag;
          lambda_iso = 0.25*exp(1.5);
          lambda_ad = 0.25*pow(3./5,1.5);
          lambda_1 = lambda_ad+(lambda_iso-lambda_ad)*pow(gamma_cooling*gamma_cooling/(88+gamma_cooling*gamma_cooling),0.22);
          lambda_2 = exp(4.5/(3+pow(beta_compton_drag,0.75)))*1/(pow(pow(1+beta_compton_drag,0.5)+1,2));
          lambda = lambda_1*lambda_2/lambda_iso;
          rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in kg/m^3 */
          M_b_dot = 4*_PI_*lambda*rho*r_B*r_B*v_eff; //in kg s^-1
          T_ion = 1.5e4*_eV_over_Kelvin_;
          tau_cooling = 1.5/(5+pow(gamma_cooling,2./3));
          Y_s = pow((1+x_e_infinity)/2,2./3*13.6/T_ion)*tau_cooling/4*pow(1-5./2*tau_cooling,1./3)*m_p/m_e;
          T_s = m_e * Y_s*pow(1+Y_s/0.27,-1./3); // in MeV
          if(T_s/m_e > 1)  J = 27/(2*_PI_)*(log(2*T_s/(m_e)*exp(-0.577)+0.08)+4./3);
          else J = 4/_PI_*sqrt(2/_PI_)*pow(T_s/m_e,-0.5)*(1+5.5*pow(T_s/m_e,1.25));
          L_ed = 4*_PI_*_G_*preco->PBH_accreting_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);
          L_acc = 1./137*T_s/(m_p)*J*pow(M_b_dot*_c_*_c_,2)/L_ed;
         }

        *energy_rate =  (rho_cdm_today/(preco->PBH_accreting_mass*M_sun*_c_*_c_))*pow(1+z,3)*L_acc*preco->PBH_fraction;
        // fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e\n",z, beta_compton_drag, gamma_cooling,lambda,M_b_dot*_c_*_c_/L_ed,T_s*1e6/_eV_over_Kelvin_,T_s*J/m_p/137,v_eff,v_B,v_l,L_acc_2/L_ed,*energy_rate);
        // fprintf(stdout, "%e %e %e %e %e %e %e %e \n",x_e, M_b_dot,lambda,m_dot_2,l2,L_acc_2,*energy_rate,z);
        // fprintf(stdout, "%e %e %e \n", z,M_b_dot*_c_*_c_/L_ed,L_acc_2/L_ed);
        free(pvecback);

}

//GFA, this function computes the energy injection coming from accreting black holes
// for different values of the PBH masses, and store the results in preco->energy_rate_at_mass[index_M]
int thermodynamics_accreting_pbh_energy_injection_PBH_MF(
                                                         struct precision * ppr,
                                                         struct background * pba,
                                                         struct recombination * preco,
                                                         double z,
                                                         ErrorMsg error_message
                                                        ){

  double rho_cdm_today;
  double tau;
  int last_index_back;
  int index_M; // GFA
  double * pvecback;
  //Parameters related to PBH
  double c_s, v_eff,v_eff_2,v_l, r_B,x_e,beta,beta_eff,beta_hat,x_cr,lambda,n_gas,M_b_dot,M_sun,M_ed_dot,epsilon,L_acc,Integrale,Normalization;
  double m_H, m_dot, m_dot_2, L_acc_2,L_ed,l,l2,M_crit;
  double rho, m_p = 938, m_e = 0.511, T_infinity = 0, rho_infinity = 0, x_e_infinity = 0, P_infinity = 0, rho_cmb = 0, t_B = 0, v_B = 0;
  double lambda_1,lambda_2,lambda_ad,lambda_iso,gamma_cooling,beta_compton_drag, T_s, T_ion, Y_s, J,tau_cooling;
  double Value_min, Value_med, Value_max, a=0, epsilon_0=0.1;
  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */

  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  for (index_M=0; index_M < preco->num_PBH_accreting_mass; index_M++) {

    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               error_message);


          c_s = 5.7e3*pow(preco->Tm_tmp/2730,0.5);//conversion km en m
          M_sun = 2e30; // in Kg
          n_gas = 200*1e6*pow((1+z)/1000,3); // 1e6 = conversion cm^-3 en m^-3;
          m_H= 1.67e-27; // Hydrogen mass in kg

          // x_e = 1;
          x_e = preco->xe_tmp;
          T_infinity = preco->Tm_tmp*_eV_over_Kelvin_*1e-6; //Temperature in MeV
          /** Disk accretion from Poulin et al. 1707.04206 */
          if(preco->PBH_accretion_recipe == disk_accretion){
              L_ed = 4*_PI_*_G_*preco->table_PBH_accreting_mass[index_M]*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);
              M_ed_dot= 10*L_ed/(_c_*_c_);
              M_crit = 0.01*M_ed_dot;
              v_B = sqrt((1+x_e)*T_infinity/m_p)*_c_;



              if(preco->PBH_relative_velocities < 0.){
                v_l = 30*MIN(1,z/1000)*1e3; // in m/s.
                if(v_B < v_l) v_eff = sqrt(v_B*v_l);
                else v_eff = v_B;
              }
              else{
                v_l = preco->PBH_relative_velocities*1e3; // converted to m/s.
                v_eff = pow(v_l*v_l+v_B*v_B,0.5);
              }

              lambda = preco->PBH_accretion_eigenvalue;
              rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in kg/m^3 */
              M_b_dot = 4*_PI_*lambda*pow(_G_*preco->table_PBH_accreting_mass[index_M]*M_sun,2)*rho*pow(v_eff,-3.);

              if(preco->PBH_ADAF_delta == 1e-3){
                Value_min = 7.6e-5;
                Value_med = 4.5e-3;
                Value_max = 7.1e-3;
                if(M_b_dot/M_ed_dot <= Value_min){
                  epsilon_0 = 0.065;
                  a = 0.71;
                }
                else if(Value_min < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot  <= Value_med){
                  epsilon_0 = 0.020;
                  a = 0.47;
                }
                else if(Value_med < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot <= Value_max){
                  epsilon_0 = 0.26;
                  a = 3.67;
                }
                else{
                  epsilon_0 = 0.1;
                  a = 0.;
                }
                epsilon = epsilon_0 * pow(M_b_dot / M_crit,a);
              }
              else if (preco->PBH_ADAF_delta == 0.1){
                Value_min = 9.4e-5;
                Value_med = 5e-3;
                Value_max = 6.6e-3;
                if(M_b_dot/M_ed_dot <= Value_min){
                  epsilon_0 = 0.12;
                  a = 0.59;
                }
                else if(Value_min < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot  <= Value_med){
                  epsilon_0 = 0.026;
                  a = 0.27;
                }
                else if(Value_med < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot <= Value_max){
                  epsilon_0 = 0.50;
                  a = 4.53;
                }
                else{
                  epsilon_0 = 0.1;
                  a = 0.;
                }
                epsilon = epsilon_0 * pow(M_b_dot / M_crit,a);
              }
              else if (preco->PBH_ADAF_delta == 0.5){

                Value_min = 2.9e-5;
                Value_med = 3.3e-3;
                Value_max = 5.3e-3;
                if(M_b_dot/M_ed_dot <= Value_min){
                  epsilon_0 = 1.58;
                  a = 0.65;
                }
                else if(Value_min < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot  <= Value_med){
                  epsilon_0 = 0.055;
                  a = 0.076;
                }
                else if(Value_med < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot <= Value_max){
                  epsilon_0 = 0.17;
                  a = 1.12;
                }
                else{
                  epsilon_0 = 0.1;
                  a = 0.;
                }
                epsilon = epsilon_0 * pow(M_b_dot / M_crit,a);
              }

              L_acc = epsilon*M_b_dot*_c_*_c_;

            }
          /** Spherical accretion from Ali-Haimoud et al. 1612.05644 */
          else if(preco->PBH_accretion_recipe == spherical_accretion){
            rho_cmb = pvecback[pba->index_bg_rho_g]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_*_c_*_c_* 6.241509e12; /* energy density in MeV/m^3 */
            // x_e_infinity = 1; // change to 1 for the strong-feedback case
            x_e_infinity = x_e; // change to x_e for the no-feedback case
            v_B = sqrt((1+x_e_infinity)*T_infinity/m_p)*_c_; //sound speed.
            if(preco->PBH_relative_velocities < 0.){
              v_l = 30*MIN(1,z/1000)*1e3; // in m/s.
              if(v_B < v_l) v_eff = sqrt(v_B*v_l);
              else v_eff = v_B;
            }
            else{
              v_l = preco->PBH_relative_velocities*1e3; // converted to m/s.
              v_eff = pow(v_l*v_l+v_B*v_B,0.5);
            }
            r_B = _G_*preco->table_PBH_accreting_mass[index_M]*M_sun*pow(v_eff,-2); // in m
            t_B = _G_*preco->table_PBH_accreting_mass[index_M]*M_sun/pow(v_eff,3); // in s
            beta_compton_drag = 4./3*x_e_infinity*_sigma_*rho_cmb*t_B/(m_p)*_c_;
            gamma_cooling = 2*m_p/(m_e*(1+x_e_infinity))*beta_compton_drag;
            lambda_iso = 0.25*exp(1.5);
            lambda_ad = 0.25*pow(3./5,1.5);
            lambda_1 = lambda_ad+(lambda_iso-lambda_ad)*pow(gamma_cooling*gamma_cooling/(88+gamma_cooling*gamma_cooling),0.22);
            lambda_2 = exp(4.5/(3+pow(beta_compton_drag,0.75)))*1/(pow(pow(1+beta_compton_drag,0.5)+1,2));
            lambda = lambda_1*lambda_2/lambda_iso;
            rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in kg/m^3 */
            M_b_dot = 4*_PI_*lambda*rho*r_B*r_B*v_eff; //in kg s^-1
            T_ion = 1.5e4*_eV_over_Kelvin_;
            tau_cooling = 1.5/(5+pow(gamma_cooling,2./3));
            Y_s = pow((1+x_e_infinity)/2,2./3*13.6/T_ion)*tau_cooling/4*pow(1-5./2*tau_cooling,1./3)*m_p/m_e;
            T_s = m_e * Y_s*pow(1+Y_s/0.27,-1./3); // in MeV
            if(T_s/m_e > 1)  J = 27/(2*_PI_)*(log(2*T_s/(m_e)*exp(-0.577)+0.08)+4./3);
            else J = 4/_PI_*sqrt(2/_PI_)*pow(T_s/m_e,-0.5)*(1+5.5*pow(T_s/m_e,1.25));
            L_ed = 4*_PI_*_G_*preco->table_PBH_accreting_mass[index_M]*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);
            L_acc = 1./137*T_s/(m_p)*J*pow(M_b_dot*_c_*_c_,2)/L_ed;
           }

            preco->energy_rate_at_mass[index_M] = (L_acc/preco->table_PBH_accreting_mass[index_M]);
         //    preco->energy_rate_at_mass[index_M] = (rho_cdm_today/(preco->table_PBH_accreting_mass[index_M]*M_sun*_c_*_c_))*pow(1+z,3)*L_acc*preco->PBH_fraction;

          if (isinf(preco->energy_rate_at_mass[index_M])==1) {
            printf("energy_rate_at_mass is inf at z =%e\n",z);
            exit(EXIT_FAILURE);
          }

  }
  free(pvecback);

  return _SUCCESS_;

}



/*************************New version, improved by Vivian Poulin******************************/
int thermodynamics_onthespot_energy_injection(
                                              struct precision * ppr,
                                              struct background * pba,
                                              struct recombination * preco,
                                              double z,
                                              double * energy_rate,
                                              ErrorMsg error_message
                                              ) {
  //fprintf(stdout,"Level: >>thermodynamics_on_the_spot_energy_injection<<< | z = %e | energy rate (before) = %e |",z,*energy_rate);
  if(preco->annihilation > 0){
    thermodynamics_DM_annihilation_energy_injection(ppr,pba,preco,z,energy_rate,error_message);
  }
  if(preco->decay_fraction > 0.){
    thermodynamics_DM_decay_energy_injection(ppr,pba,preco,z,energy_rate,error_message);
  }
  if(preco->PBH_accreting_mass > 0.){
    thermodynamics_accreting_pbh_energy_injection(ppr,pba,preco,z,energy_rate,error_message);
  }
  if(preco->PBH_evaporating_mass > 0.){
    thermodynamics_evaporating_pbh_energy_injection(ppr,pba,preco,z,energy_rate,error_message);
  }
  //fprintf(stdout," energy rate (after) = %e \n",*energy_rate);
  /* energy density rate in J/m^3/s (remember that annihilation_at_z is in m^3/s/Kg and decay in s^-1) */
  return _SUCCESS_;

}

/**
 * In case of non-minimal cosmology, this function determines the
 * effective energy rate absorbed by the IGM at a given redshift
 * (beyond the on-the-spot annihilation). This energy injection may
 * come e.g. from dark matter annihilation or decay.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param preco Input: pointer to recombination structure
 * @param z Input: redshift
 * @param energy_rate Output: energy density injection rate
 * @param error_message Output: error message
 * @return the error status
 */

int thermodynamics_energy_injection(
                                    struct precision * ppr,
                                    struct background * pba,
                                    struct recombination * preco,
                                    double z,
                                    double * energy_rate,
                                    ErrorMsg error_message
                                    ) {

  double zp,dz;
  double integrand,first_integrand;
  double factor,result;
  double nH0;
  double onthespot;
  double exponent_z,exponent_zp;
  if (preco->annihilation > 0 || preco->decay_fraction > 0 || preco->PBH_accreting_mass > 0 || preco->PBH_evaporating_mass > 0 ) {
    //fprintf(stdout,"Level >>thermodynamics_energy_injection<< | z = %e | energy rate = %e\n",z,*energy_rate);
    if (preco->has_on_the_spot == _FALSE_) {


      if(preco->energy_deposition_function == Analytical_approximation){
        // // /* number of hydrogen nuclei today in m**-3 */
        nH0 = 3.*preco->H0*preco->H0*pba->Omega0_b/(8.*_PI_*_G_*_m_H_)*(1.-preco->YHe);

        /*Value from Poulin et al. 1508.01370*/
        /* factor = c sigma_T n_H(0) / (H(0) \sqrt(Omega_m)) (dimensionless) */
        // factor = _sigma_ * nH0 / pba->H0 * _Mpc_over_m_ / sqrt(pba->Omega0_b+pba->Omega0_cdm);
        // exponent_z = 8;
        // exponent_zp = 7.5;

        /*Value from Ali-Haimoud & Kamionkowski 1612.05644*/
        factor = 0.1*_sigma_ * nH0 / pba->H0 * _Mpc_over_m_ / sqrt(pba->Omega0_b+pba->Omega0_cdm);
        exponent_z = 7;
        exponent_zp = 6.5;


        /* integral over z'(=zp) with step dz */
        dz=1.;

        /* first point in trapezoidal integral */
        zp = z;
        class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,zp,&onthespot,error_message),
                   error_message,
                   error_message);
        first_integrand = factor*pow(1+z,exponent_z)/pow(1+zp,exponent_zp)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 7 and 7.5
        result = 0.5*dz*first_integrand;

        /* other points in trapezoidal integral */
        do {

          zp += dz;
          class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,zp,&onthespot,error_message),
                     error_message,
                     error_message);
          integrand = factor*pow(1+z,exponent_z)/pow(1+zp,exponent_zp)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 7 and 7.5
          result += dz*integrand;

        } while (integrand/first_integrand > 0.02);
        if(result < 1e-100) result=0.;
      }

      // // /***********************************************************************************************************************/
      else if(preco->energy_deposition_function == function_from_file){

            if(preco->energy_repart_coefficient!=no_factorization){
              class_call(thermodynamics_annihilation_f_eff_interpolate(ppr,pba,preco,z,error_message),
                        error_message,
                        error_message);
              preco->f_eff=MAX(preco->f_eff,0.);
            }
            else preco->f_eff=1.;

            class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,z,&result,error_message),
                      error_message,
                      error_message);
            result =  result*preco->f_eff;
            // fprintf(stdout, "energy_rate %e preco->f_eff %e\n", result,preco->f_eff);
      }
      // // /***********************************************************************************************************************/
      else if(preco->energy_deposition_function == DarkAges){

            class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,z,&result,error_message),
                      error_message,
                      error_message);
            // fprintf(stdout, "energy_rate %e \n", result);
      }

      // /* uncomment these lines if you also want to compute the on-the-spot for comparison */
      // class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,z,&onthespot,error_message),
      //            error_message,
      //            error_message);
      // /* these test lines print the energy rate rescaled by (1+z)^6 in J/m^3/s, with or without the on-the-spot approximation */
      // fprintf(stdout,"%e  %e  %e  %e\n",
      // 1.+z,
      // result/pow(1.+z,6),
      // onthespot/pow(1.+z,6),result/onthespot);
    }
    else {
      class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,z,&result,error_message),
                 error_message,
                 error_message);
      if(preco->f_eff>0)result *= preco->f_eff; //If preco->f_eff is defined, here we multiply by f_eff.
      // fprintf(stdout, "energy_rate %e preco->f_eff %e\n", result,preco->f_eff);

       /* effective energy density rate in J/m^3/s  */
    }
    *energy_rate = result;

  }
  else {
    *energy_rate = 0.;
  }

  return _SUCCESS_;

}

/**
 * This subroutine contains the reionization function \f$ X_e(z) \f$
 * (one for each scheme; so far, only the function corresponding to
 * the reio_camb scheme is coded)
 *
 * @param z     Input: redshift
 * @param pth   Input: pointer to thermo structure, to know which scheme is used
 * @param preio Input: pointer to reionization structure, containing the parameters of the function \f$ X_e(z) \f$
 * @param xe    Output: \f$ X_e(z) \f$
 */

int thermodynamics_reionization_function(
                                         double z,
                                         struct thermo * pth,
                                         struct reionization * preio,
                                         struct recombination * preco,
                                         double * xe
                                         ) {
  /** Summary: */
  /** - define local variables */
  double argument, A, B, factor;
  int i;
  double z_jump;
  int jump;
  double center,before, after,width,one_jump;
  double x_tmp,z_tmp;
  /** - implementation of ionization function similar to the one in CAMB */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
    /** -> case z > z_reio_start */

    if(z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];

    }

    else {

      /** - --> case z < z_reio_start: hydrogen contribution (tanh of complicated argument) */
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */

      A = (pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                    preio->reionization_parameters[preio->index_reio_exponent])
                - pow((1.+z),preio->reionization_parameters[preio->index_reio_exponent]));
      // A = MAX (A,0.);
      B = (preio->reionization_parameters[preio->index_reio_exponent]
        *pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
             (preio->reionization_parameters[preio->index_reio_exponent]-1.)))*preio->reionization_parameters[preio->index_reio_width];
             if(B == 0)argument = 0;
             else argument = A / B;



      if (pth->reio_parametrization == reio_camb) {
        *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])
          *(tanh(argument)+1.)/2.
          +preio->reionization_parameters[preio->index_reio_xe_before];
      }

      else {
        *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])
          *tanh(argument)
          +preio->reionization_parameters[preio->index_reio_xe_before];
      }

        /** -> case z < z_reio_start: helium contribution (tanh of simpler argument) */

          argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
            /preio->reionization_parameters[preio->index_helium_fullreio_width];
          /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
          *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
            * (tanh(argument)+1.)/2.;

      if(*xe > 0.1*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_10_percent == 0){
        pth->z_10_percent = z;
        // fprintf(stdout, "pth->z_10_percent %e\n", pth->z_10_percent);
      }
      if(*xe > 0.50*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_50_percent == 0 ){
        pth->z_50_percent = z;
        // fprintf(stdout, "pth->z_50_percent %e\n", pth->z_50_percent);

      }
      if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent == 0){
      // if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent > 200){
        pth->z_99_percent = z;
        // class_test(pth->z_99_percent<=6,
        //            pth->error_message,
        //            "z_99_percent < 6, we reject the point");

      }
      if(pth->z_99_percent != 0 && pth->z_10_percent != 0 && pth->duration_of_reionization == 0){
        pth->duration_of_reionization = pth->z_10_percent  - pth->z_99_percent;
        if(pth->duration_of_reionization < 0) pth->duration_of_reionization ==0;
        // class_test(pth->duration_of_reionization<1,
        //            pth->error_message,
        //            "duration_of_reionization < 1, we reject the point");
        // fprintf(stdout, "pth->duration_of_reionization %e v2 %e \n", pth->duration_of_reionization,pth->z_10_percent  -preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]);

      }
      x_tmp = *xe;
      z_tmp = z;
    }

    return _SUCCESS_;

  }

  if(pth->reio_parametrization == reio_douspis_et_al){
    if(z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];

    }

    else {

    if(z < preio->reionization_parameters[preio->index_zp_douspis_et_al]){
      factor = (1-preio->reionization_parameters[preio->index_Qp_douspis_et_al])/(pow(1+preio->reionization_parameters[preio->index_zp_douspis_et_al],3)-1)
      *(pow(1+preio->reionization_parameters[preio->index_zp_douspis_et_al],3)-pow(1+z,3))
      +preio->reionization_parameters[preio->index_Qp_douspis_et_al];

    *xe =  (preio->reionization_parameters[preio->index_reio_xe_after]
             -preio->reionization_parameters[preio->index_reio_xe_before])*
             factor;

    }
    else{
      factor = preio->reionization_parameters[preio->index_Qp_douspis_et_al]*exp(-preio->reionization_parameters[preio->index_lambda_douspis_et_al]
        *pow(z-preio->reionization_parameters[preio->index_zp_douspis_et_al],3)/(pow(z-preio->reionization_parameters[preio->index_zp_douspis_et_al],2)+0.2));
      *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])*
               factor;
      *xe = MAX(preio->reionization_parameters[preio->index_reio_xe_before],*xe);

    }

    /** -> case z < z_reio_start: helium contribution (tanh of simpler argument) */
    /********Helium reionization is taken into account with a camb-like parametrization***********/
      argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
        /preio->reionization_parameters[preio->index_helium_fullreio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
      *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
        * (tanh(argument)+1.)/2.;
    // fprintf(stdout, "z %e x_e %e xe_before %e factor %elambda %e zp %e qp %e\n", z,*xe,preio->reionization_parameters[preio->index_reio_xe_before],factor,preio->reionization_parameters[preio->index_lambda_douspis_et_al],preio->reionization_parameters[preio->index_zp_douspis_et_al],preio->reionization_parameters[preio->index_Qp_douspis_et_al]);

    }
    if(*xe > 0.1*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_10_percent == 0){
      pth->z_10_percent = z;
      // fprintf(stdout, "pth->z_10_percent %e\n", pth->z_10_percent);
    }
    if(*xe > 0.50*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_50_percent == 0 ){
      pth->z_50_percent = z;
      // fprintf(stdout, "pth->z_50_percent %e\n", pth->z_50_percent);

    }
    if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent == 0){
    // if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent > 200){
      pth->z_99_percent = z;
    //   class_test(pth->z_99_percent<=6,
    //              pth->error_message,
    //              "z_99_percent < 6, we reject the point");
    //  fprintf(stdout, "pth->z_99_percent %e  \n", pth->z_99_percent);


    }
    if(pth->z_99_percent != 0 && pth->z_10_percent != 0 && pth->duration_of_reionization == 0){
      pth->duration_of_reionization = pth->z_10_percent  - pth->z_99_percent;
      if(pth->duration_of_reionization < 0) pth->duration_of_reionization ==0;
      // class_test(pth->duration_of_reionization<1,
      //            pth->error_message,
      //            "duration_of_reionization < 1, we reject the point");
      // fprintf(stdout, "pth->duration_of_reionization %e v2 %e \n", pth->duration_of_reionization,pth->z_10_percent  -preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]);

    }

    return _SUCCESS_;

  }
  if(pth->reio_parametrization == reio_asymmetric_planck_16){
    if(z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];

    }

    else {

    if(z < preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]){


    *xe =  (preio->reionization_parameters[preio->index_reio_xe_after]
             -preio->reionization_parameters[preio->index_reio_xe_before]);

    }
    else{
      factor = pow((preio->reionization_parameters[preio->index_reio_start]-z)/(preio->reionization_parameters[preio->index_reio_start]-preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]),preio->reionization_parameters[preio->index_alpha_asymmetric_planck_16]);
      *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])*
               factor;
      *xe = MAX(preio->reionization_parameters[preio->index_reio_xe_before],*xe);

    }

    if(*xe > 0.1*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_10_percent == 0){
      pth->z_10_percent = z;
      // fprintf(stdout, "pth->z_10_percent %e\n", pth->z_10_percent);
    }
    if(*xe > 0.50*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_50_percent == 0 ){
      pth->z_50_percent = z;
      // fprintf(stdout, "pth->z_50_percent %e\n", pth->z_50_percent);

    }
    if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent == 0){
    // if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent > 200){
      pth->z_99_percent = z;
      // class_test(pth->z_99_percent<=6,
      //            pth->error_message,
      //            "z_99_percent < 6, we reject the point");


    }
    if(pth->z_99_percent != 0 && pth->z_10_percent != 0 && pth->duration_of_reionization == 0){
      pth->duration_of_reionization = pth->z_10_percent  - pth->z_99_percent;
      if(pth->duration_of_reionization < 0) pth->duration_of_reionization ==0;
      // class_test(pth->duration_of_reionization<1,
      //            pth->error_message,
      //            "duration_of_reionization < 1, we reject the point");
      // fprintf(stdout, "pth->duration_of_reionization %e v2 %e \n", pth->duration_of_reionization,pth->z_10_percent  -preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]);

    }
    // fprintf(stdout, "z %e xe %e factor %e\n",z,*xe,factor);

    /** -> case z < z_reio_start: helium contribution (tanh of simpler argument) */
    /********Helium reionization is taken into account with a camb-like parametrization***********/
      argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
        /preio->reionization_parameters[preio->index_helium_fullreio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
      *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
        * (tanh(argument)+1.)/2.;
    // fprintf(stdout, "z %e x_e %e xe_before %e factor %elambda %e zp %e qp %e\n", z,*xe,preio->reionization_parameters[preio->index_reio_xe_before],factor,preio->reionization_parameters[preio->index_lambda_douspis_et_al],preio->reionization_parameters[preio->index_zp_douspis_et_al],preio->reionization_parameters[preio->index_Qp_douspis_et_al]);



    }


    return _SUCCESS_;

  }
  if(pth->reio_parametrization == reio_stars_sfr_source_term){
    if(z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];
    }

    else {

    /** -> case z < z_reio_start: helium contribution (tanh of simpler argument) */
    /********Helium reionization is taken into account with a camb-like parametrization***********/
      argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
        /preio->reionization_parameters[preio->index_helium_fullreio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
      *xe =  (preio->reionization_parameters[preio->index_reio_xe_after]
             -preio->reionization_parameters[preio->index_reio_xe_before])
        *(tanh(argument)+1.)/2.
        +preio->reionization_parameters[preio->index_reio_xe_before]; //both He reionization are done simulatenouesly
      // *xe= MIN(*xe,preio->reionization_parameters[preio->index_reio_xe_after]);
      // fprintf(stderr, "z %e xe %e  xe after %e xe before %e Yhe %e argument %e \n",z, *xe,preio->reionization_parameters[preio->index_reio_xe_after],preio->reionization_parameters[preio->index_reio_xe_before],preio->reionization_parameters[preio->index_helium_fullreio_fraction],(tanh(argument)+1.)/2);
    }
    return _SUCCESS_;

  }

  /** - implementation of binned ionization function similar to astro-ph/0606552 */

  if (pth->reio_parametrization == reio_bins_tanh) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1];
    }

    else if (z < preio->reionization_parameters[preio->index_reio_first_z]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe];
    }

    else {

      i = 0;
      while (preio->reionization_parameters[preio->index_reio_first_z+i+1]<z) i++;
           /* This is the expression of the tanh-like jumps of the
              reio_bins_tanh scheme until the 10.06.2015. It appeared to be
              not robust enough. It could lead to a kink in xe(z) near the
              maximum value of z at which reionisation is sampled. It has
              been replaced by the simpler and more robust expression
              below.

             *xe = preio->reionization_parameters[preio->index_reio_first_xe+i]
               +0.5*(tanh((2.*(z-preio->reionization_parameters[preio->index_reio_first_z+i])
                           /(preio->reionization_parameters[preio->index_reio_first_z+i+1]
      @@ -1412,27 +1907,7 @@ int thermodynamics_reionization_function(
                     /tanh(1./preio->reionization_parameters[preio->index_reio_step_sharpness])+1.)
               *(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
                 -preio->reionization_parameters[preio->index_reio_first_xe+i]);
           */

           /* compute the central redshift value of the tanh jump */

           if (i == preio->reio_num_z-2) {
             z_jump = preio->reionization_parameters[preio->index_reio_first_z+i]
               + 0.5*(preio->reionization_parameters[preio->index_reio_first_z+i]
                      -preio->reionization_parameters[preio->index_reio_first_z+i-1]);
           }
           else  {
             z_jump =  0.5*(preio->reionization_parameters[preio->index_reio_first_z+i+1]
                            + preio->reionization_parameters[preio->index_reio_first_z+i]);
           }

           /* implementation of the tanh jump */

           *xe = preio->reionization_parameters[preio->index_reio_first_xe+i]
             +0.5*(tanh((z-z_jump)
                        /preio->reionization_parameters[preio->index_reio_step_sharpness])+1.)
             *(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
               -preio->reionization_parameters[preio->index_reio_first_xe+i]);


    }

    return _SUCCESS_;

  }

  /** - implementation of many tanh jumps */

  if (pth->reio_parametrization == reio_many_tanh) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1];
    }

    else if (z > preio->reionization_parameters[preio->index_reio_first_z]) {

      *xe = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1];

      for (jump=1; jump<preio->reio_num_z-1; jump++){

        center = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1-jump];
        // before and after are meant with respect to growing z, not growing time
        before = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1-jump]
          -preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-jump];
        after = 0.;
        width = preio->reionization_parameters[preio->index_reio_step_sharpness];

        class_call(thermodynamics_tanh(z,center,before,after,width,&one_jump),
                   pth->error_message,
                   pth->error_message);

        *xe += one_jump;

      }

    }

    else {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe];
    }

    return _SUCCESS_;

  }

    /** - implementation of reio_inter */

  if (pth->reio_parametrization == reio_inter) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1];
      class_stop(pth->error_message,"Check: is it normal that we are here?");
    }

    else {

      i=0;
      while (preio->reionization_parameters[preio->index_reio_first_z+i+1] < z) i++;

      double z_min = preio->reionization_parameters[preio->index_reio_first_z+i];
      double z_max = preio->reionization_parameters[preio->index_reio_first_z+i+1];

      class_test(z<z_min,
                 pth->error_message,
                 "");

      class_test(z>z_max,
                 pth->error_message,
                 "");

      double x=(z-preio->reionization_parameters[preio->index_reio_first_z+i])
        /(preio->reionization_parameters[preio->index_reio_first_z+i+1]
          -preio->reionization_parameters[preio->index_reio_first_z+i]);

      *xe = preio->reionization_parameters[preio->index_reio_first_xe+i]
        + x*(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
             -preio->reionization_parameters[preio->index_reio_first_xe+i]);

      class_test(*xe<0.,
                 pth->error_message,
                 "%e %e %e\n",
                 x,
                 preio->reionization_parameters[preio->index_reio_first_xe+i],
                 preio->reionization_parameters[preio->index_reio_first_xe+i+1]);

    }

    return _SUCCESS_;

  }
  class_test(0 == 0,
             pth->error_message,
             "value of reio_parametrization=%d unclear",pth->reio_parametrization);
}

/**
 * This subroutine reads \f$ X_e(z) \f$ in the recombination table at
 * the time at which reionization starts. Hence it provides correct
 * initial conditions for the reionization function.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pth   Input: pointer to thermo structure
 * @param preco Input: pointer to recombination structure
 * @param z     Input: redshift z_reio_start
 * @param xe    Output: \f$ X_e(z) \f$ at z
 */

int thermodynamics_get_xe_before_reionization(
                                              struct precision * ppr,
                                              struct thermo * pth,
                                              struct recombination * preco,
                                              double z,
                                              double * xe
                                              ) {

  int last_index=0;

  class_call(array_interpolate_one_growing_closeby(preco->recombination_table,
                                                   preco->re_size,
                                                   preco->rt_size,
                                                   preco->index_re_z,
                                                   z,
                                                   &last_index,
                                                   preco->index_re_xe,
                                                   xe,
                                                   pth->error_message),
             pth->error_message,
             pth->error_message);

  return _SUCCESS_;

}


/**
 * This routine computes the reionization history. In the reio_camb
 * scheme, this is straightforward if the input parameter is the
 * reionization redshift. If the input is the optical depth, need to
 * find z_reio by dichotomy (trying several z_reio until the correct
 * tau_reio is approached).
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermo structure
 * @param preco Input: pointer to filled recombination structure
 * @param preio Input/Output: pointer to reionization structure (to be filled)
 * @param pvecback   Input: vector of background quantities (used as workspace: must be already allocated, with format short_info or larger, but does not need to be filled)
 * @return the error status
 */

int thermodynamics_reionization(
                                struct precision * ppr,
                                struct background * pba,
                                struct thermo * pth,
                                struct recombination * preco,
                                struct reionization * preio,
                                double * pvecback
                                ) {

  /** Summary: */

  /** - define local variables */

  int counter;
  double z_sup,z_mid,z_inf;
  double tau_sup,tau_mid,tau_inf;
  int bin;
  int point;
  double xe_input,xe_actual;

  /** - allocate the vector of parameters defining the function \f$ X_e(z) \f$ */



  class_alloc(preio->reionization_parameters,preio->reio_num_params*sizeof(double),pth->error_message);

  /** (a) if reionization implemented like in CAMB */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh) ) {

    /** - --> set values of these parameters, excepted those depending on the reionization redshift */

    if (pth->reio_parametrization == reio_camb ) {
      preio->reionization_parameters[preio->index_reio_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */
    }
    if (pth->reio_parametrization == reio_half_tanh ) {
      preio->reionization_parameters[preio->index_reio_xe_after] = 1.;
      // ; /* xe_after_reio: neglect He ionization */
      // + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized H */
      // + 2*pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + fully ionized He */
    }
    preio->reionization_parameters[preio->index_reio_exponent] = pth->reionization_exponent; /* reio_exponent */
    preio->reionization_parameters[preio->index_reio_width] = pth->reionization_width;    /* reio_width */
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */

    class_test(preio->reionization_parameters[preio->index_reio_exponent]==0,
               pth->error_message,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_reio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");


    /** - --> if reionization redshift given as an input, initialize the remaining values and fill reionization table*/

    if (pth->reio_z_or_tau == reio_z) {

      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = pth->z_reio;

      /* infer starting redshift for hydrogen */

      if (pth->reio_parametrization == reio_camb || pth->reio_parametrization == reio_half_tanh) {

        preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*pth->reionization_width;

        /* if starting redshift for helium is larger, take that one
           (does not happen in realistic models) */
        if (preio->reionization_parameters[preio->index_reio_start] <
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width)

          preio->reionization_parameters[preio->index_reio_start] =
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;

      }
      else {

        preio->reionization_parameters[preio->index_reio_start] = pth->z_reio;
      }

      class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                 pth->error_message,
                 "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

      /* infer xe_before_reio */
      class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                           pth,
                                                           preco,
                                                           preio->reionization_parameters[preio->index_reio_start],
                                                           &(preio->reionization_parameters[preio->index_reio_xe_before])),
                 pth->error_message,
                 pth->error_message);


      /* fill reionization table */
      class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                 pth->error_message,
                 pth->error_message);

      pth->tau_reio=preio->reionization_optical_depth;

    }

    /** - --> if reionization optical depth given as an input, find reionization redshift by dichotomy and initialize the remaining values */

    if (pth->reio_z_or_tau == reio_tau) {


      /* upper value */

      z_sup = ppr->reionization_z_start_max-ppr->reionization_start_factor*pth->reionization_width;
      class_test(z_sup < 0.,
                 pth->error_message,
                 "parameters are such that reionization cannot take place before today while starting after z_start_max; need to increase z_start_max");

      /* maximum possible reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_sup;
      /* maximum possible starting redshift */
      preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
      /* infer xe_before_reio */

      class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                           pth,
                                                           preco,
                                                           preio->reionization_parameters[preio->index_reio_start],
                                                           &(preio->reionization_parameters[preio->index_reio_xe_before])),
                 pth->error_message,
                 pth->error_message);

      /* fill reionization table */
      class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                 pth->error_message,
                 pth->error_message);


      tau_sup=preio->reionization_optical_depth;

      class_test(tau_sup < pth->tau_reio,
                 pth->error_message,
                 "parameters are such that reionization cannot start after z_start_max");

      /* lower value */

      z_inf = 0.;
      tau_inf = 0.;

      /* try intermediate values */

      counter=0;
      while ((tau_sup-tau_inf) > pth->tau_reio * ppr->reionization_optical_depth_tol) {
        z_mid=0.5*(z_sup+z_inf);

        /* reionization redshift */
        preio->reionization_parameters[preio->index_reio_redshift] = z_mid;
        /* infer starting redshift for hygrogen */
        preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*pth->reionization_width;
        /* if starting redshift for helium is larger, take that one
           (does not happen in realistic models) */
        if (preio->reionization_parameters[preio->index_reio_start] <
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width)

          preio->reionization_parameters[preio->index_reio_start] =
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;

        class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                   pth->error_message,
                   "starting redshift for reionization > reionization_z_start_max = %e",ppr->reionization_z_start_max);

        /* infer xe_before_reio */
        class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                             pth,
                                                             preco,
                                                             preio->reionization_parameters[preio->index_reio_start],
                                                             &(preio->reionization_parameters[preio->index_reio_xe_before])),
                   pth->error_message,
                   pth->error_message);

        /* clean and fill reionization table */
        free(preio->reionization_table);
        class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                   pth->error_message,
                   pth->error_message);

        tau_mid=preio->reionization_optical_depth;

        /* trial */

        if (tau_mid > pth->tau_reio) {
          z_sup=z_mid;
          tau_sup=tau_mid;
        }
        else {
          z_inf=z_mid;
          tau_inf=tau_mid;
        }

        counter++;
        class_test(counter > _MAX_IT_,
                   pth->error_message,
                   "while searching for reionization_optical_depth, maximum number of iterations exceeded");
      }

      /* store z_reionization in thermodynamics structure */
      pth->z_reio=preio->reionization_parameters[preio->index_reio_redshift];

    }

    free(preio->reionization_parameters);

    return _SUCCESS_;

  }

  /** - (b) if reionization implemented with reio_bins_tanh scheme */

  if (pth->reio_parametrization == reio_bins_tanh) {

    /* this algorithm requires at least two bin centers (i.e. at least
       4 values in the (z,xe) array, counting the edges). */
    class_test(pth->binned_reio_num<2,
               pth->error_message,
               "current implementation of binned reio requires at least two bin centers");

    /* check that this input can be interpreted by the code */
    for (bin=1; bin<pth->binned_reio_num; bin++) {
      class_test(pth->binned_reio_z[bin-1]>=pth->binned_reio_z[bin],
                 pth->error_message,
                 "value of reionization bin centers z_i expected to be passed in growing order: %e, %e",
                 pth->binned_reio_z[bin-1],
                 pth->binned_reio_z[bin]);
    }

    /* the code will not only copy here the "bin centers" passed in
       input. It will add an initial and final value for (z,xe).
       First, fill all entries except the first and the last */

    for (bin=1; bin<preio->reio_num_z-1; bin++) {
      preio->reionization_parameters[preio->index_reio_first_z+bin] = pth->binned_reio_z[bin-1];
      preio->reionization_parameters[preio->index_reio_first_xe+bin] = pth->binned_reio_xe[bin-1];
    }


    /* find largest value of z in the array. We choose to define it as
       z_(i_max) + 2*(the distance between z_(i_max) and z_(i_max-1)). E.g. if
       the bins are in 10,12,14, the largest z will be 18. */
    preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1] =

      preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-2]
      +2.*(preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-2]
        -preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-3]);

    /* copy this value in reio_start */
    preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
               pth->error_message,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

    /* find smallest value of z in the array. We choose
       to define it as z_0 - (the distance between z_1 and z_0). E.g. if
       the bins are in 10,12,14, the stop redshift will be 8. */

    preio->reionization_parameters[preio->index_reio_first_z] =
      2.*preio->reionization_parameters[preio->index_reio_first_z+1]
      -preio->reionization_parameters[preio->index_reio_first_z+2];

    /* check it's not too small */
    /* 6.06.2015: changed this test to simply imposing that the first z is at least zero */
    /*
    class_test(preio->reionization_parameters[preio->index_reio_first_z] < 0,
               pth->error_message,
               "final redshift for reionization = %e, you must change the binning or redefine the way in which the code extrapolates below the first value of z_i",preio->reionization_parameters[preio->index_reio_first_z]);
    */
    if (preio->reionization_parameters[preio->index_reio_first_z] < 0) {
      preio->reionization_parameters[preio->index_reio_first_z] = 0.;
    }

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                         pth,
                                                         preco,
                                                         preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
                                                         &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
               pth->error_message,
               pth->error_message);

    /* infer xe after reio */
    preio->reionization_parameters[preio->index_reio_first_xe] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */

    /* pass step sharpness parameter */
    preio->reionization_parameters[preio->index_reio_step_sharpness] = pth->binned_reio_step_sharpness;

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
               pth->error_message,
               pth->error_message);

    pth->tau_reio=preio->reionization_optical_depth;

    return _SUCCESS_;

  }

  /** - (c) if reionization implemented with reio_many_tanh scheme */

  if (pth->reio_parametrization == reio_many_tanh) {

    /* this algorithm requires at least one jump centers */
    class_test(pth->many_tanh_num<1,
               pth->error_message,
               "current implementation of reio_many_tanh requires at least one jump center");

    /* check that z input can be interpreted by the code */
    for (bin=1; bin<pth->many_tanh_num; bin++) {
      class_test(pth->many_tanh_z[bin-1]>=pth->many_tanh_z[bin],
                 pth->error_message,
                 "value of reionization bin centers z_i expected to be passed in growing order: %e, %e",
                 pth->many_tanh_z[bin-1],
                 pth->many_tanh_z[bin]);
    }

    /* the code will not only copy here the "jump centers" passed in
       input. It will add an initial and final value for (z,xe).
       First, fill all entries except the first and the last */

    for (bin=1; bin<preio->reio_num_z-1; bin++) {

      preio->reionization_parameters[preio->index_reio_first_z+bin] = pth->many_tanh_z[bin-1];

      /* check that xe input can be interpreted by the code */
      xe_input = pth->many_tanh_xe[bin-1];
      if (xe_input >= 0.) {
        xe_actual = xe_input;
      }
      //-1 means "after hydrogen + first helium recombination"
      else if ((xe_input<-0.9) && (xe_input>-1.1)) {
        xe_actual = 1. + pth->YHe/(_not4_*(1.-pth->YHe));
      }
      //-2 means "after hydrogen + second helium recombination"
      else if ((xe_input<-1.9) && (xe_input>-2.1)) {
        xe_actual = 1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe));
      }
      //other negative number is nonsense
      else {
        class_stop(pth->error_message,
                   "Your entry for many_tanh_xe[%d] is %e, this makes no sense (either positive or 0,-1,-2)",
                   bin-1,pth->many_tanh_xe[bin-1]);
      }

      preio->reionization_parameters[preio->index_reio_first_xe+bin] = xe_actual;
    }

    /* find largest value of z in the array. We choose to define it as
       z_(i_max) + ppr->reionization_start_factor*step_sharpness. */
    preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1] =
      preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-2]
      +ppr->reionization_start_factor*pth->many_tanh_width;

    /* copy this value in reio_start */
    preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
               pth->error_message,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

    /* find smallest value of z in the array. We choose
       to define it as z_0 - ppr->reionization_start_factor*step_sharpness, but at least zero. */

    preio->reionization_parameters[preio->index_reio_first_z] =
      preio->reionization_parameters[preio->index_reio_first_z+1]
      -ppr->reionization_start_factor*pth->many_tanh_width;

    if (preio->reionization_parameters[preio->index_reio_first_z] < 0) {
      preio->reionization_parameters[preio->index_reio_first_z] = 0.;
    }

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                         pth,
                                                         preco,
                                                         preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
                                                         &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
               pth->error_message,
               pth->error_message);

    /* infer xe after reio */

    preio->reionization_parameters[preio->index_reio_first_xe] = preio->reionization_parameters[preio->index_reio_first_xe+1];

    /* if we want to model only hydrogen reionization and neglect both helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1.;

    /* if we want to model only hydrogen + first helium reionization and neglect second helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));

    /* if we want to model hydrogen + two helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe));

    /* pass step sharpness parameter */
    class_test(pth->many_tanh_width<=0,
               pth->error_message,
               "many_tanh_width must be strictly positive, you passed %e",
               pth->many_tanh_width);

    preio->reionization_parameters[preio->index_reio_step_sharpness] = pth->many_tanh_width;

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
               pth->error_message,
               pth->error_message);

    pth->tau_reio=preio->reionization_optical_depth;

    return _SUCCESS_;

  }
  if((pth->reio_parametrization == reio_douspis_et_al)){
    preio->reionization_parameters[preio->index_reio_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */
    preio->reionization_parameters[preio->index_lambda_douspis_et_al] = pth->lambda_douspis_et_al;
    preio->reionization_parameters[preio->index_zp_douspis_et_al] = pth->zp_douspis_et_al;
    preio->reionization_parameters[preio->index_Qp_douspis_et_al] = pth->Qp_douspis_et_al;

    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");
     preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;

     class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                pth->error_message,
                "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

     /* infer xe_before_reio */
     class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                          pth,
                                                          preco,
                                                          preio->reionization_parameters[preio->index_reio_start],
                                                          &(preio->reionization_parameters[preio->index_reio_xe_before])),
                pth->error_message,
                pth->error_message);
     /* fill reionization table */
     class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                pth->error_message,
                pth->error_message);
     pth->tau_reio=preio->reionization_optical_depth;

     return _SUCCESS_;

  }

  if((pth->reio_parametrization == reio_asymmetric_planck_16)){
    preio->reionization_parameters[preio->index_reio_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */
    preio->reionization_parameters[preio->index_alpha_asymmetric_planck_16] = pth->alpha_asymmetric_planck_16;
    preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16] = pth->z_end_asymmetric_planck_16;
    preio->reionization_parameters[preio->index_reio_start] = pth->z_start_asymmetric_planck_16;
    pth->z_10_percent = 0;
    pth->z_50_percent = 0;
    pth->z_99_percent = 0;
    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");
    //  preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;

     class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                pth->error_message,
                "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

     /* infer xe_before_reio */
     class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                          pth,
                                                          preco,
                                                          preio->reionization_parameters[preio->index_reio_start],
                                                          &(preio->reionization_parameters[preio->index_reio_xe_before])),
                pth->error_message,
                pth->error_message);
     /* fill reionization table */
     class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                pth->error_message,
                pth->error_message);
     pth->tau_reio=preio->reionization_optical_depth;

     return _SUCCESS_;

  }

  if((pth->reio_parametrization == reio_stars_sfr_source_term)){//Helium reionization is still a tanh, to be improved
    preio->reionization_parameters[preio->index_reio_xe_after] = 1. + 2*pth->YHe/(_not4_*(1.-pth->YHe));
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */
    preio->reionization_parameters[preio->index_reio_start] = pth->helium_fullreio_redshift+pth->helium_fullreio_width;
    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");
    //  preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
     class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                pth->error_message,
                "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

     /* infer xe_before_reio */
     class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                          pth,
                                                          preco,
                                                          preio->reionization_parameters[preio->index_reio_start],
                                                          &(preio->reionization_parameters[preio->index_reio_xe_before])),
                pth->error_message,
                pth->error_message);
     /* fill reionization table */
     class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                pth->error_message,
                pth->error_message);
     pth->tau_reio=preio->reionization_optical_depth;

     return _SUCCESS_;

  }
  /** - (d) if reionization implemented with reio_inter scheme */

  if (pth->reio_parametrization == reio_inter) {

    /* this parametrization requires at least one point (z,xe) */
    class_test(pth->reio_inter_num<1,
               pth->error_message,
               "current implementation of reio_inter requires at least one point (z,xe)");

    /* this parametrization requires that the first z value is zero */
    class_test(pth->reio_inter_z[0] != 0.,
               pth->error_message,
               "For reio_inter scheme, the first value of reio_inter_z[...]  should always be zero, you passed %e",
               pth->reio_inter_z[0]);

    /* check that z input can be interpreted by the code */
    for (point=1; point<pth->reio_inter_num; point++) {
      class_test(pth->reio_inter_z[point-1]>=pth->reio_inter_z[point],
                 pth->error_message,
                 "value of reionization bin centers z_i expected to be passed in growing order, unlike: %e, %e",
                 pth->reio_inter_z[point-1],
                 pth->reio_inter_z[point]);
    }

    /* this parametrization requires that the last x_i value is zero
       (the code will substitute it with the value that one would get in
       absence of reionization, as compute by the recombination code) */
    class_test(pth->reio_inter_xe[pth->reio_inter_num-1] != 0.,
               pth->error_message,
               "For reio_inter scheme, the last value of reio_inter_xe[...]  should always be zero, you passed %e",
               pth->reio_inter_xe[pth->reio_inter_num-1]);

    /* copy here the (z,xe) values passed in input. */

    for (point=0; point<preio->reio_num_z; point++) {

      preio->reionization_parameters[preio->index_reio_first_z+point] = pth->reio_inter_z[point];

      /* check that xe input can be interpreted by the code */
      xe_input = pth->reio_inter_xe[point];
      if (xe_input >= 0.) {
        xe_actual = xe_input;
      }
      //-1 means "after hydrogen + first helium recombination"
      else if ((xe_input<-0.9) && (xe_input>-1.1)) {
        xe_actual = 1. + pth->YHe/(_not4_*(1.-pth->YHe));
      }
      //-2 means "after hydrogen + second helium recombination"
      else if ((xe_input<-1.9) && (xe_input>-2.1)) {
        xe_actual = 1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe));
      }
      //other negative number is nonsense
      else {
        class_stop(pth->error_message,
                   "Your entry for reio_inter_xe[%d] is %e, this makes no sense (either positive or 0,-1,-2)",
                   point,pth->reio_inter_xe[point]);
      }

      preio->reionization_parameters[preio->index_reio_first_xe+point] = xe_actual;
    }

    /* copy highest redshift in reio_start */
    preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
               pth->error_message,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                         pth,
                                                         preco,
                                                         preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
                                                         &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
               pth->error_message,
               pth->error_message);

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
               pth->error_message,
               pth->error_message);

    pth->tau_reio=preio->reionization_optical_depth;

    return _SUCCESS_;

  }
  class_test(0 == 0,
             pth->error_message,
             "value of reio_z_or_tau=%d unclear",pth->reio_z_or_tau);
}

/**
 * For fixed input reionization parameters, this routine computes the
 * reionization history and fills the reionization table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermo structure
 * @param preco Input: pointer to filled recombination structure
 * @param preio Input/Output: pointer to reionization structure (to be filled)
 * @param pvecback   Input: vector of background quantities (used as workspace: must be already allocated, with format short_info or larger, but does not need to be filled)
 * @return the error status
 */

int thermodynamics_reionization_sample(
                                       struct precision * ppr,
                                       struct background * pba,
                                       struct thermo * pth,
                                       struct recombination * preco,
                                       struct reionization * preio,
                                       double * pvecback
                                       ) {

  /** Summary: */

  /** - define local variables */

  /* a growing table (since the number of redshift steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* current vector of values related to reionization */
  double * reio_vector;
  /* running index inside thermodynamics table */
  int i,j;
  int number_of_redshifts;
  // GFA
  int k, index_M;
  double energy_rate_dep_heat;

  /* values of z, dz, X_e */
  double dz,dz_max;
  double z,z_next;
  double xe,xe_next,x_tmp;
  double dkappadz,dkappadz_next;
  double delta_z_old, delta_z_new;
  double Tb,Yp,dTdz,dTdz_adia,dTdz_CMB,dTdz_DM,dTdz_stars,opacity,mu;
  double dkappadtau,dkappadtau_next;
  double energy_rate;
  double tau;
  double chi_heat, chi_heat_x_ray;
  double chi_lya;
  double chi_ionH;
  double chi_ionHe;
  double chi_lowE;
  double argument;
  int last_index_back;
  double relative_variation;
  double L_x, rho_sfr;
  double rho_cdm_today; // GFA
  double M_sun = 2e30; // in Kg

  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */


  Yp = pth->YHe;

  /** - (a) allocate vector of values related to reionization */
  class_alloc(reio_vector,preio->re_size*sizeof(double),pth->error_message);

  /** - (b) create a growTable with gt_init() */
  class_call(gt_init(&gTable),
             gTable.error_message,
             pth->error_message);



  /** - (c) first line is taken from thermodynamics table, just before reionization starts */

  /** - --> look where to start in current thermodynamics table */
  i=0;
  while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < preio->reionization_parameters[preio->index_reio_start]) {
    i++;
    class_test(i == ppr->recfast_Nz0,
               pth->error_message,
               "reionization_z_start_max = %e > largest redshift in thermodynamics table",ppr->reionization_z_start_max);
  }
  if(preco->recombination_table[i*preco->re_size+preco->index_re_z] >  preio->reionization_parameters[preio->index_reio_start])i--;
  j=i;
  /** - --> get redshift */
  z=preco->recombination_table[i*preco->re_size+preco->index_re_z];

  reio_vector[preio->index_re_z]=z;
  preio->index_reco_when_reio_start=i;

  /** - --> get \f$ X_e \f$ */
  class_call(thermodynamics_reionization_function(z,pth,preio,preco,&xe),
             pth->error_message,
             pth->error_message);

    if(pth->reio_stars_and_dark_matter == _TRUE_){
      xe=preco->recombination_table[i*preco->re_size+preco->index_re_xe];
      // xe=MAX(xe,x_tmp);
    }
  reio_vector[preio->index_re_xe] = xe;

  /** -  --> get \f$ d \kappa / d z = (d \kappa / d \tau) * (d \tau / d z) = - (d \kappa / d \tau) / H \f$ */

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pth->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->short_info,
                               pba->inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             pth->error_message);

  reio_vector[preio->index_re_dkappadtau] = (1.+z) * (1.+z) * pth->n_e * xe * _sigma_ * _Mpc_over_m_;

  class_test(pvecback[pba->index_bg_H] == 0.,
             pth->error_message,
             "stop to avoid division by zero");

  reio_vector[preio->index_re_dkappadz] = reio_vector[preio->index_re_dkappadtau] / pvecback[pba->index_bg_H];

  dkappadz = reio_vector[preio->index_re_dkappadz];
  dkappadtau = reio_vector[preio->index_re_dkappadtau];


  /** - --> get baryon temperature **/
  Tb = preco->recombination_table[i*preco->re_size+preco->index_re_Tb];
  reio_vector[preio->index_re_Tb] = Tb;

  /** - --> after recombination, Tb scales like (1+z)**2. Compute constant factor Tb/(1+z)**2. */
  //Tba2 = Tb/(1+z)/(1+z);

  /** - --> get baryon sound speed */
  reio_vector[preio->index_re_cb2] = 5./3. * _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * Yp + xe * (1.-Yp)) * Tb;

  /** - --> store these values in growing table */
  class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
             gTable.error_message,
             pth->error_message);

  number_of_redshifts=1;

  /** - (d) set the maximum step value (equal to the step in thermodynamics table) */
  dz_max=preco->recombination_table[i*preco->re_size+preco->index_re_z]
    -preco->recombination_table[(i-1)*preco->re_size+preco->index_re_z];
    // fprintf(stderr, "dz %e\n", dz_max);
  /** - (e) loop over redshift values in order to find values of z, x_e, kappa' (Tb and cb2 found later by integration). The sampling in z space is found here. */

  /* initial step */
  dz = dz_max;
  while (z > 0.) {
    if (j<0)j=0;
    // fprintf(stdout, "j %d \n",j);
    // dz = MAX(ppr->smallest_allowed_variation,dz);

    class_test(dz < ppr->smallest_allowed_variation,
               pth->error_message,
               "stuck in the loop for reionization sampling, as if you were trying to impose a discontinuous evolution for xe(z)");

    /* - try next step */
    z_next=z-dz;
    delta_z_old = z_next-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z];
    delta_z_new = z_next-preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z];
    if(fabs(delta_z_old)<fabs(delta_z_new))j++;
    while(z_next > preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z])j++;
    // fprintf(stdout, "z = %e z_next = %e\n",z,z_next);

    if (z_next < 0.) z_next=0.;
    class_call(thermodynamics_reionization_function(z_next,pth,preio,preco,&xe_next),
               pth->error_message,
               pth->error_message);

    if(pth->reio_stars_and_dark_matter == _TRUE_){
      /**
       * This small routine compares the reionization table to the recombination one and choose the highest x_e between the two.
       * This way allows to enables to avoid unphysical discontiniuty in the ionization fraction at low x_e.
       * First, we interpolate the value of x_e at the evaluated redshift from the recombination table. The linear interpolation has been checked
       * to work well but it could be improve for security. Then we perform the comparison.
       */

      x_tmp= (preco->recombination_table[(j-2)*preco->re_size+preco->index_re_xe]-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_xe])/(preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z]
        -preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z])*(z_next-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z])+
        preco->recombination_table[(j-1)*preco->re_size+preco->index_re_xe]  ;
      x_tmp = MAX(0.,x_tmp); // Small check to avoid negative values of x_e.

      if(x_tmp <1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe))) xe_next=MAX(xe_next,x_tmp); // Here the comparison is made.
      else x_tmp = 1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe)); // the maximal value that x_e can reach.

    }

    class_call(background_tau_of_z(pba,
                                   z_next,
                                   &tau),
               pba->error_message,
               pth->error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->short_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

    class_test(pvecback[pba->index_bg_H] == 0.,
               pth->error_message,
               "stop to avoid division by zero");
    // if(xe_next > 1.17) fprintf(stdout, "error xe next %e\n", xe_next);
    dkappadz_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_ / pvecback[pba->index_bg_H];
    dkappadtau_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_;


    class_test((dkappadz == 0.) || (dkappadtau == 0.),
               pth->error_message,
               "stop to avoid division by zero");

    relative_variation = fabs((dkappadz_next-dkappadz)/dkappadz) +
      fabs((dkappadtau_next-dkappadtau)/dkappadtau);

    if (relative_variation < ppr->reionization_sampling || pth->reio_stars_and_dark_matter == _TRUE_) {
      /* accept the step: get \f$ z, X_e, d kappa / d z \f$ and store in growing table */

      z=z_next;
      xe=xe_next;
      dkappadz=dkappadz_next;
      dkappadtau= dkappadtau_next;
      class_test((dkappadz == 0.) || (dkappadtau == 0.),
                 pth->error_message,
                 "dkappadz=%e, dkappadtau=%e, stop to avoid division by zero",dkappadz,dkappadtau);

      reio_vector[preio->index_re_z] = z;
      reio_vector[preio->index_re_xe] = xe;
      reio_vector[preio->index_re_dkappadz] = dkappadz;
      reio_vector[preio->index_re_dkappadtau] = dkappadz * pvecback[pba->index_bg_H];

      class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
                 gTable.error_message,
                 pth->error_message);

      number_of_redshifts++;
      // delta_z_old = z_next-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z];
      // delta_z_new = z_next-preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z];
      // if(fabs(delta_z_old)>fabs(delta_z_new))j--;
      j--;
      dz = MIN(0.9*(ppr->reionization_sampling/relative_variation),5.)*dz;
      // dz = MIN(dz,dz_max);
      // dz = MAX(ppr->smallest_allowed_variation,dz);
    }
    else {
      /* do not accept the step and update dz */
      // delta_z_old = z_next-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z];
      // delta_z_new = z_next-preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z];
      dz = 0.9*(ppr->reionization_sampling/relative_variation)*dz;
      // dz = MIN(dz,dz_max);
      // dz = MAX(ppr->smallest_allowed_variation,dz);
      // j--;
      // if(fabs(delta_z_old)>fabs(delta_z_new))j--;

    }

  }

  /** - (f) allocate reionization_table with correct size */
  class_alloc(preio->reionization_table,preio->re_size*number_of_redshifts*sizeof(double),pth->error_message);

  preio->rt_size=number_of_redshifts;

  /** - (g) retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             pth->error_message);

  /** - (h) copy growTable to reionization_temporary_table (invert order of lines, so that redshift is growing, like in recombination table) */
  for (i=0; i < preio->rt_size; i++) {
    memcopy_result = memcpy(preio->reionization_table+i*preio->re_size,pData+(preio->rt_size-i-1)*preio->re_size,preio->re_size*sizeof(double));
    class_test(memcopy_result != preio->reionization_table+i*preio->re_size,
               pth->error_message,
               "cannot copy data back to reionization_temporary_table");

  }

  /** - (i) free the growTable with gt_free() , free vector of reionization variables */
  class_call(gt_free(&gTable),
             gTable.error_message,
             pth->error_message);

  free(reio_vector);

  /** - (j) another loop on z, to integrate equation for Tb and to compute cb2 */
  for (i=preio->rt_size-1; i >0 ; i--) {
    z = preio->reionization_table[i*preio->re_size+preio->index_re_z];
    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               pth->error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

    dz = (preio->reionization_table[i*preio->re_size+preio->index_re_z]-preio->reionization_table[(i-1)*preio->re_size+preio->index_re_z]);

    opacity = (1.+z) * (1.+z) * pth->n_e
      * preio->reionization_table[i*preio->re_size+preio->index_re_xe] * _sigma_ * _Mpc_over_m_;

    mu = _m_H_/(1. + (1./_not4_ - 1.) * pth->YHe + preio->reionization_table[i*preio->re_size+preio->index_re_xe] * (1.-pth->YHe));


    /** - derivative of baryon temperature */

      /** - First possibility: Add a tanh term in the temperature and bypass evolution equation*/
      if(pth->star_heating_parametrization== heating_reiolike_tanh){

        argument = (pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                        preio->reionization_parameters[preio->index_reio_exponent])
                    - pow((1.+z),preio->reionization_parameters[preio->index_reio_exponent]))
          /(preio->reionization_parameters[preio->index_reio_exponent]
            /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
            *pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                 (preio->reionization_parameters[preio->index_reio_exponent]-1.)))
          /preio->reionization_parameters[preio->index_reio_width];
        /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
          if(z< preio->reionization_parameters[preio->index_reio_start])
          preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb] = (pth->final_IGM_temperature
                 -preio->reionization_table[i*preio->re_size+preio->index_re_Tb])
        *(tanh(argument)+1.)/2 + preio->reionization_table[i*preio->re_size+preio->index_re_Tb];

      }
      /** - Second possibility: Compute temperature evolution from each sources*/

      else {


    dTdz_adia=2./(1+z)*preio->reionization_table[i*preio->re_size+preio->index_re_Tb];

    dTdz_CMB = - 2.*mu/_m_e_*4.*pvecback[pba->index_bg_rho_g]/3./pvecback[pba->index_bg_rho_b]*opacity*
      (pba->T_cmb * (1.+z)-preio->reionization_table[i*preio->re_size+preio->index_re_Tb])/pvecback[pba->index_bg_H];

      /** - Parameters related to exotic energy injection */
      if((pth->annihilation != 0 || pth->decay_fraction != 0 || pth->PBH_accreting_mass != 0 || pth->PBH_evaporating_mass != 0 || pth->has_extended_PBH_MassFunc == _TRUE_)){ // GFA

            /** - --> derivative of baryon temperature */
              preco->xe_tmp=preio->reionization_table[i*preio->re_size+preio->index_re_xe];
              preco->Tm_tmp=preio->reionization_table[i*preio->re_size+preio->index_re_Tb];

              if (pth->has_extended_PBH_MassFunc == _TRUE_) { // GFA, compute energy injection per each pbh mass
                class_call(thermodynamics_accreting_pbh_energy_injection_PBH_MF(ppr,pba,preco,z,pth->error_message),
                           pth->error_message,
                           pth->error_message);
              } else {
                class_call(thermodynamics_energy_injection(ppr,pba,preco,z,&energy_rate,pth->error_message),
                           pth->error_message,
                           pth->error_message);
              }


              preco->z_tmp=z;
               /* coefficient as revised by Galli et al. 2013 (in fact it is an interpolation by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013) */
               if(pth->energy_repart_coefficient==GSVI || pth->energy_repart_coefficient ==chi_from_file){
                 class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,preio->reionization_table[i*preio->re_size+preio->index_re_xe]),
                          pth->error_message,
                          pth->error_message);
                 chi_heat = pth->chi_heat;
              }
               if(pth->energy_repart_coefficient==no_factorization){
                 if (pth->has_extended_PBH_MassFunc == _TRUE_) { // GFA, interpolate f(z) per channel and now also per each mass
                   class_call(thermodynamics_annihilation_coefficients_interpolate_PBH_MF(ppr,pba,pth,z),
                            pth->error_message,
                            pth->error_message);
                 } else {
                   class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,z),
                            pth->error_message,
                            pth->error_message);
                   chi_heat = pth->chi_heat;
                 }

              }
               /* old approximation from Chen and Kamionkowski */
               if(pth->energy_repart_coefficient==SSCK){
                 chi_heat = (1.+2.*preio->reionization_table[i*preio->re_size+preio->index_re_xe])/3.;
               }

              chi_heat= MIN(chi_heat,1.);
              chi_heat = MAX(chi_heat,0.);
              if (pth->has_extended_PBH_MassFunc == _TRUE_) {
              for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) {
              pth->chi_heat_at_mass[index_M]= MIN(pth->chi_heat_at_mass[index_M],1.);
              pth->chi_heat_at_mass[index_M] = MAX(pth->chi_heat_at_mass[index_M],0.);
              }
              }

              if (pth->has_extended_PBH_MassFunc == _TRUE_) {
               for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) { // GFA, evaluate integrand at each mass
                 pth->ener_rate_dep_heat_per_mass[index_M] = pth->table_PBH_MassFunc[index_M]*pth->chi_heat_at_mass[index_M]*preco->energy_rate_at_mass[index_M];
               }
               energy_rate_dep_heat = 0.;
               for (k=0; k < pth->num_PBH_accreting_mass-1; k++) { //integral over the mass using trapezoidal rule, improve to Simpson, although with at least 50 mass bins it shouldn't matter
                 energy_rate_dep_heat += 0.5*(pth->table_PBH_accreting_mass[k+1]-pth->table_PBH_accreting_mass[k])*(pth->ener_rate_dep_heat_per_mass[k]+pth->ener_rate_dep_heat_per_mass[k+1]);
               }
                energy_rate_dep_heat *= preco->PBH_fraction*(rho_cdm_today/(M_sun*_c_*_c_))*pow(1+z,3);
             } else { //standard case, no integral required
                energy_rate_dep_heat = energy_rate*chi_heat;
               }

              dTdz_DM = - 2./(3.*_k_B_)*energy_rate_dep_heat
              /(preco->Nnow*pow(1.+z,3))/(1.+preco->fHe+preio->reionization_table[i*preio->re_size+preio->index_re_xe])
              /(pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_*(1.+z)); /* energy injection */

              if(pth->thermodynamics_verbose>10){
                fprintf(stdout, "z %e dTdz_CMB %e dTdz_DM %e xe %e energy_rate %e chi_heat %e\n",z,dTdz_CMB,dTdz_DM,preio->reionization_table[i*preio->re_size+preio->index_re_xe],energy_rate, chi_heat);
              }

      }
      else dTdz_DM = 0.;


      /** parametrization of reheating by stars */
      if(pth->star_heating_parametrization == heating_none){ //Standard assumption, no reheating by stars. Enough for accurate computations of CMB power spectra.
        dTdz_stars = 0;
      }

      else if(pth->star_heating_parametrization== heating_stars_sfr_source_term ){ //Reheating term based on the SFR rates. See Poulin et al. 1508.01370.
          rho_sfr = pth->ap*pow(1+z,pth->bp)/(1+pow((1+z)/pth->cp,pth->dp))/pow(_Mpc_over_m_,3)*pow(1+z,3)*(1+tanh((pth->z_start_reio_stars-z)))/2; //add a (sharp) smoothing function.
          L_x = pth->Ex* pth->fx *rho_sfr*2./(3.*_k_B_*preco->Nnow*pow(1.+z,3)*(1.+preco->fHe+preio->reionization_table[i*preio->re_size+preio->index_re_xe]))
          /(pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_*(1.+z));

          dTdz_stars = -L_x*(1+2*preio->reionization_table[i*preio->re_size+preio->index_re_xe])/3.;
        }

      // else if(){
      //   // ready for other parametrization
      // }

      dTdz = dTdz_adia+dTdz_CMB+dTdz_DM+dTdz_stars;
      if(pth->thermodynamics_verbose>10){
      // fprintf(stdout, "z %e dT %e Tmat %e dTdz_adia %e dTdz_CMB %e dTdz_DM %e dTdz_stars %e opacity %e xe %e energy_rate %e chi_heat %e pvecback[pba->index_bg_H] %e\n", z,dTdz, preio->reionization_table[i*preio->re_size+preio->index_re_Tb],dTdz_adia, dTdz_CMB ,dTdz_DM,energy_rate,chi_heat,pvecback[pba->index_bg_H]);
      fprintf(stdout, "z %e dT %e Tmat %e dTdz_adia %e dTdz_CMB %e dTdz_DM %e dTdz_stars %e opacity %e xe %e energy_rate %e chi_heat %e H %e\n", z,dTdz, preio->reionization_table[i*preio->re_size+preio->index_re_Tb],dTdz_adia, dTdz_CMB ,dTdz_DM,dTdz_stars,opacity,preio->reionization_table[i*preio->re_size+preio->index_re_xe],energy_rate,chi_heat,pvecback[pba->index_bg_H]);
      }
      /** - --> increment baryon temperature  */

        preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb] =
          preio->reionization_table[i*preio->re_size+preio->index_re_Tb]-dTdz*dz;

      }


    /** - get baryon sound speed */

    preio->reionization_table[(i-1)*preio->re_size+preio->index_re_cb2] = _k_B_/ ( _c_ * _c_ * mu)
      * preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb]
      *(1.+(1+z)/3.*dTdz/preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb]);

  }


  /** - --> spline \f$ d \tau / dz \f$ with respect to z in view of integrating for optical depth */
  class_call(array_spline(preio->reionization_table,
                          preio->re_size,
                          preio->rt_size,
                          preio->index_re_z,
                          preio->index_re_dkappadz,
                          preio->index_re_d3kappadz3,
                          _SPLINE_EST_DERIV_,
                          pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> integrate for optical depth */
  class_call(array_integrate_all_spline(preio->reionization_table,
                                        preio->re_size,
                                        preio->rt_size,
                                        preio->index_re_z,
                                        preio->index_re_dkappadz,
                                        preio->index_re_d3kappadz3,
                                        &(preio->reionization_optical_depth),
                                        pth->error_message),
             pth->error_message,
             pth->error_message);


  return _SUCCESS_;

}

/**
 * Integrate thermodynamics with your favorite recombination code.
 *
 */

int thermodynamics_recombination(
                                 struct precision * ppr,
                                 struct background * pba,
                                 struct thermo * pth,
                                 struct recombination * preco,
                                 double * pvecback
                                 ) {

  if (pth->recombination==hyrec) {

    class_call(thermodynamics_recombination_with_hyrec(ppr,pba,pth,preco,pvecback),
               pth->error_message,
               pth->error_message);

  }

  if (pth->recombination==recfast) {

    class_call(thermodynamics_recombination_with_recfast(ppr,pba,pth,preco,pvecback),
               pth->error_message,
               pth->error_message);

  }

  if (pth->recombination==cosmorec) {
        class_call(thermodynamics_recombination_with_cosmorec(ppr,pba,pth,preco,pvecback),
               pth->error_message,
               pth->error_message);
  }

  return _SUCCESS_;

}
/**
 * Integrate thermodynamics with CosmoRec.
 *
 * Integrate thermodynamics with CosmoRec, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_recombination(), from thermodynamics_init().
 *
 *************************************************************************************************
 *                 CosmoRec: Cosmological Recombination Project
 *                Written by Jens Chluba (University of Manchester)
 *************************************************************************************************
 * @param ppr      Input: pointer to precision structure
 * @param pba      Input: pointer to background structure
 * @param pth      Input: pointer to thermodynamics structure
 * @param preco    Output: pointer to recombination structure
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 */
int thermodynamics_recombination_with_cosmorec(
                                            struct precision * ppr,
                                            struct background * pba,
                                            struct thermo * pth,
                                            struct recombination * preco,
                                            double * pvecback
                                            ) {
#ifdef COSMOREC
  int i;
  double nH0 = 11.223846333047*pba->Omega0_b*pba->h*pba->h*(1.-pth->YHe);  /* number density of hydrogen today in m-3 */
  double rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*pba->Omega0_cdm*_c_*_c_; /* energy density in J/m^3 */
  double DM_annihilation =  pth->annihilation*1e-6/_c_/_c_*pow(rho_cdm_today,2)/nH0*1e6/_eV_; /*conversion in cosmorec unit as described in Chluba 2010 0910.3663 (without factor 2, to respect class convention of majorana particles)*/
  preco->f_eff = 1;

  if(preco->energy_deposition_function == function_from_file){

        if(preco->energy_repart_coefficient!=no_factorization){
          class_call(thermodynamics_annihilation_f_eff_interpolate(ppr,pba,preco,600),
                    pth->error_message,
                    pth->error_message);
          preco->f_eff=MAX(preco->f_eff,0.);
        }
        else{
          class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,600),
                         pth->error_message,
                         pth->error_message);
          preco->f_eff = (pth->chi_heat+pth->chi_ionH+pth->chi_ionHe+pth->chi_lya); // we use the corrected scheme which will multiply the SSCK prescription (currently hardcoded in cosmorec).
        }
  }
  // // /***********************************************************************************************************************/
  else if(preco->energy_deposition_function == DarkAges){
      class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,600),
                     pth->error_message,
                     pth->error_message);
      preco->f_eff = (pth->chi_heat+pth->chi_ionH+pth->chi_ionHe+pth->chi_lya); // we use the corrected scheme which will multiply chen&kamionkowski's prescription (currently hardcoded in cosmorec).
  }



  DM_annihilation *= preco->f_eff;
  double runpars[4] = {
    DM_annihilation, /* defines the dark matter annihilation efficiency in eV/s. */
    pth->cosmorec_accuracy, /* setting for cosmorec accuracy (default = default cosmorec setting) */
    pth->cosmorec_verbose, /* setting for cosmorec verbose (default = no output produced) */
    pth->Lambda_over_theoritical_Lambda *_Lambda_ /* theoritical value by Labzowsky et al 2005 for H1_A2s_1s is rescaled, by default Lambda_over_theoritical_Lambda = 1. In agreement with standard cosmorec.*/
  };

  double H0 = pba->H0 / 1e3 * _c_;
  int nz = ppr->recfast_Nz0;
  double * z_arr;
  double * Hz_arr;
  double z, xe, Tm, Hz;
  double z_start=ppr->recfast_z_initial;
  double z_end=0;
  double step;

  double tau_at_z;
  int last_index;

  double * xe_out;
  double * tb_out;

  double drho_dt = 0, Tg;
  double dlnTb_dz;
  int label=0; /* iterator for cosmorec output file name, not used in class version of cosmorec */

  /* Initialize Hubble rate for CosmoRec */
  class_alloc(z_arr, sizeof(double) * nz, pth->error_message);
  class_alloc(Hz_arr, sizeof(double) * nz, pth->error_message);

  step = (z_start - z_end) / (nz);
  for(i=0; i < nz; i++) {
    z_arr[i] = z_end + i * step;

      class_call(
        background_tau_of_z(
          pba,
          z_arr[i],
          &tau_at_z
        ),
        pba->error_message,
        pth->error_message
      );

      class_call(
        background_at_tau(
          pba,
          tau_at_z,
          pba->short_info,
          pba->inter_normal,
          &last_index,
          pvecback
        ),
        pba->error_message,
        pth->error_message
      );

      Hz_arr[i]=pvecback[pba->index_bg_H] * _c_ / _Mpc_over_m_;
  }
  /* Initialize x_e and tb output tables */

  class_alloc(xe_out, sizeof(double) * nz, pth->error_message);
  class_alloc(tb_out, sizeof(double) * nz, pth->error_message);

  /* call cosmorec */
  /* Currently we give parameters separetely, eventually to be changed for a structure, easier to modify.*/

  cosmorec_calc_h_cpp_(
    &(pth->cosmorec_runmode), runpars,
    &(pba->Omega0_cdm), &(pba->Omega0_b), &(pba->Omega0_k),
    &(pba->Neff), &H0,
    &(pba->T_cmb), &(pth->YHe),
    z_arr, Hz_arr, &nz,
    z_arr, xe_out, tb_out,
    &nz,
    &label
  );

  /** - fill a few parameters in preco and pth */



  preco->rt_size = nz;
  preco->H0 = pba->H0 * _c_ / _Mpc_over_m_;
  /* preco->H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
  preco->YHe = pth->YHe;
  preco->Nnow = 3.*preco->H0*preco->H0*pba->Omega0_b*(1.-preco->YHe)/(8.*_PI_*_G_*_m_H_);
  /* energy injection parameters */
  preco->annihilation = pth->annihilation;
  preco->has_on_the_spot = pth->has_on_the_spot;
  preco->decay_fraction = pth->decay_fraction;
  preco->annihilation_f_halo = pth->annihilation_f_halo;
  preco->annihilation_z_halo = pth->annihilation_z_halo;
  pth->n_e=preco->Nnow;

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) and fill it */

  class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);

  for(i=nz-1; i >= 0; i--) {

    /** - --> get redshift, corresponding results from cosmorec, and background quantities */

    z = z_arr[i];
    xe = xe_out[i];
    Tm = tb_out[i];
    Hz = Hz_arr[i];
    /** - --> store the results in the table */

    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* redshift */
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_z)=z;

    /* ionization fraction */
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_xe)=xe;

    /* Tb */
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_Tb)=Tm;

    /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)
       with (1+z)dlnTb/dz= - [dlnTb/dlna] */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */

    Tg = pba->T_cmb * (1+z);
    dlnTb_dz = - Tg/Tm*drho_dt/(1+z)/Hz+1/(1+z);

   evaluate_TM(z, xe,preco->fHe, Tm/Tg, Tg, Hz, &drho_dt);
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_cb2)
      = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + xe * (1.-pth->YHe)) * Tm * (1. + (1+z)*dlnTb_dz / 3.);
    /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_dkappadtau)
      = (1.+z) * (1.+z) * preco->Nnow * xe * _sigma_ * _Mpc_over_m_;
      //  fprintf(stdout,"xe %e Tm %e cb2 %e z %e dlnTb_dz %e *dkappa_dtau %e\n",xe,Tm,*(preco->recombination_table+(i)*preco->re_size+preco->index_re_cb2),z,dlnTb_dz,*(preco->recombination_table+(i)*preco->re_size+preco->index_re_dkappadtau));

  }

  /* clean up */

  free(xe_out);
  free(tb_out);

  free(z_arr);
  free(Hz_arr);

#else

class_stop(pth->error_message,
           "you compiled without including the CosmoRec code, and now wish to use it. Either set the input parameter 'recombination' to something else than 'CosmoRec', or recompile after setting in the Makefile the appropriate path COSMOREC=... ");


#endif /* COSMOREC */

  return _SUCCESS_;

}

/**
 * Integrate thermodynamics with HyRec.
 *
 * Integrate thermodynamics with HyRec, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_recombination(), from thermodynamics_init().
 *
 *************************************************************************************************
 *                 HYREC: Hydrogen and Helium Recombination Code
 *         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)
 *************************************************************************************************
 *
 *
 * @param ppr      Input: pointer to precision structure
 * @param pba      Input: pointer to background structure
 * @param pth      Input: pointer to thermodynamics structure
 * @param preco    Output: pointer to recombination structure
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 */

 int thermodynamics_recombination_with_hyrec(
                                             struct precision * ppr,
                                             struct background * pba,
                                             struct thermo * pth,
                                             struct recombination * preco,
                                             double * pvecback
                                             ) {
   /** Summary: */
 #ifdef HYREC

   HYREC_DATA hyrec_data;
   hyrec_allocate(&hyrec_data, ppr->recfast_z_initial, 0.);

   double Omega_m = pba->Omega0_b + pba->Omega0_cdm + pba->Omega0_ncdm_tot;

   double alpha_ratio = 1.;    /* Ratio of fine-structure constant to standard value */
   double me_ratio    = 1.;    /* Ratio of electron mass to standard value */

   double pann        = 1.78266e-21 *pth->annihilation;  /* Converting from m^3/s/kg to cm^3/s/GeV */
   double pann_halo   = 1.78266e-21 *pth->annihilation_f_halo;

   int i,j,Nz;
   double z, xe, Tm, Hz;
   void * buffer;
   double tau;
   int last_index_back,last_index_back2;
   int on_the_spot = 1;

   if(pth->has_on_the_spot == _FALSE_){
     on_the_spot = 0;
   }


   /** - Fill the recombination structure with all important parameters */
   class_call(fill_recombination_structure(ppr,pba,pth,preco),pth->error_message,pth->error_message);

   /** - Compute the recombination history by calling hyrec_compute.
         No CLASS-like error management here, but YAH working on it :) **/

           hyrec_data.cosmo->h = pba->h;
           hyrec_data.cosmo->T0 = pba->T_cmb;
           hyrec_data.cosmo->orh2  = 4.48162687719e-7 *pba->T_cmb*pba->T_cmb*pba->T_cmb*pba->T_cmb *(1. + 0.227107317660239 * pba->Neff);
           hyrec_data.cosmo->obh2 = pba->Omega0_b*pba->h*pba->h;
           hyrec_data.cosmo->omh2 = Omega_m*pba->h*pba->h;
           hyrec_data.cosmo->inj_params->odmh2 = hyrec_data.cosmo->omh2 - hyrec_data.cosmo->obh2;
           hyrec_data.cosmo->okh2 = pba->Omega0_k*pba->h*pba->h;
           hyrec_data.cosmo->odeh2 = (1.-Omega_m - hyrec_data.cosmo->okh2 - hyrec_data.cosmo->orh2/pba->h/pba->h)*pba->h*pba->h;
          //  hyrec_data.w0 = pba->w0_fld;
          //  hyrec_data.wa = pba->wa_fld;
           hyrec_data.cosmo->Y = pth->YHe;
           hyrec_data.cosmo->Nnueff = pba->Neff;
           hyrec_data.cosmo->nH0 = 11.223846333047e-6*hyrec_data.cosmo->obh2*(1.-hyrec_data.cosmo->Y);  /* number density of hydrogen today in m-3 */
           hyrec_data.cosmo->fHe = hyrec_data.cosmo->Y/(1-hyrec_data.cosmo->Y)/3.97153;              /* abundance of helium by number */

           hyrec_data.cosmo->inj_params->Omega0_b = pba->Omega0_b;
           hyrec_data.cosmo->inj_params->Omega0_cdm = pba->Omega0_cdm;
           hyrec_data.cosmo->inj_params->Omega0_r = hyrec_data.cosmo->orh2/pow(pba->h,2);
           hyrec_data.cosmo->inj_params->H0 = pba->h*100;
          //  hyrec_data.cosmo->inj_params->zstart = ppr->recfast_z_initial; /* Redshift range */
          //  hyrec_data.cosmo->inj_params->zend = 0.;
          //  hyrec_data.cosmo->inj_params->dlna = 8.49e-5;
          //  hyrec_data.cosmo->inj_params->nz = (long) floor(2+log((1.+hyrec_data.cosmo->inj_params->zstart)/(1.+hyrec_data.cosmo->inj_params->zend))/hyrec_data.cosmo->inj_params->dlna);
           hyrec_data.cosmo->inj_params->pann = 1.78266e-21*pth->annihilation;
           hyrec_data.cosmo->inj_params->pann_halo = 1.78266e-21*pth->annihilation_f_halo;
           hyrec_data.cosmo->inj_params->on_the_spot = on_the_spot;
           hyrec_data.cosmo->inj_params->decay_fraction = pth->decay_fraction;
           hyrec_data.cosmo->inj_params->Gamma_dcdm = pba->Gamma_dcdm;
           hyrec_data.cosmo->inj_params->f_eff = pth->f_eff;
           hyrec_data.cosmo->inj_params->ann_f_halo = pth->annihilation_f_halo;
           hyrec_data.cosmo->inj_params->ann_z_halo = pth->annihilation_z_halo;
           // GFA
           hyrec_data.cosmo->inj_params->has_UCMH_spike = pth->has_UCMH_spike;
           hyrec_data.cosmo->inj_params->Number_z =ppr->Number_z;
           hyrec_data.cosmo->inj_params->z_table_for_boost = pth->z_table_for_boost;
           hyrec_data.cosmo->inj_params->boost_table = pth->boost_table;
           hyrec_data.cosmo->inj_params->annihil_coef_num_lines = pth->annihil_coef_num_lines;
           hyrec_data.cosmo->inj_params->annihil_coef_xe = pth->annihil_coef_xe;
           hyrec_data.cosmo->inj_params->annihil_coef_heat = pth->annihil_coef_heat;
           hyrec_data.cosmo->inj_params->annihil_coef_ionH = pth->annihil_coef_ionH;
           hyrec_data.cosmo->inj_params->annihil_coef_ionHe = pth->annihil_coef_ionHe;
           hyrec_data.cosmo->inj_params->annihil_coef_lya = pth->annihil_coef_lya;
           hyrec_data.cosmo->inj_params->annihil_coef_lowE = pth->annihil_coef_lowE;

           hyrec_data.cosmo->inj_params->annihil_coef_dd_heat = pth->annihil_coef_dd_heat;
           hyrec_data.cosmo->inj_params->annihil_coef_dd_ionH = pth->annihil_coef_dd_ionH;
           hyrec_data.cosmo->inj_params->annihil_coef_dd_ionHe = pth->annihil_coef_dd_ionHe;
           hyrec_data.cosmo->inj_params->annihil_coef_dd_lya = pth->annihil_coef_dd_lya;
           hyrec_data.cosmo->inj_params->annihil_coef_dd_lowE = pth->annihil_coef_lowE;
           hyrec_data.cosmo->inj_params->annihil_f_eff_num_lines = preco->annihil_f_eff_num_lines;
           hyrec_data.cosmo->inj_params->annihil_z = preco->annihil_z;
           hyrec_data.cosmo->inj_params->annihil_f_eff = preco->annihil_f_eff;
           hyrec_data.cosmo->inj_params->annihil_dd_f_eff = preco->annihil_dd_f_eff;

           if ((pth->PBH_table_is_initialized) == _FALSE_ && pth->PBH_evaporating_mass > 0.) {
             pth->PBH_table_is_initialized = _TRUE_;
             preco->PBH_table_is_initialized = _TRUE_;
             PBH_evaporating_mass_time_evolution(ppr,pba,preco,pth->error_message);
             hyrec_data.cosmo->inj_params->PBH_table_is_initialized= preco->PBH_table_is_initialized;
             hyrec_data.cosmo->inj_params->PBH_table_z = preco->PBH_table_z;
             hyrec_data.cosmo->inj_params->PBH_table_mass = preco->PBH_table_mass;
             hyrec_data.cosmo->inj_params->PBH_table_mass_dd = preco->PBH_table_mass_dd;
             hyrec_data.cosmo->inj_params->PBH_table_F = preco->PBH_table_F;
             hyrec_data.cosmo->inj_params->PBH_table_F_dd = preco->PBH_table_F_dd;
             hyrec_data.cosmo->inj_params->PBH_table_size= preco->PBH_table_size;
           }

            hyrec_data.cosmo->inj_params->PBH_low_mass = pth->PBH_evaporating_mass;


           if(pth->energy_deposition_function==Analytical_approximation) {
             hyrec_data.cosmo->inj_params->energy_deposition_treatment = 0;
           }
           else if(pth->energy_deposition_function==function_from_file || pth->energy_deposition_function==DarkAges) {
             hyrec_data.cosmo->inj_params->energy_deposition_treatment = 1;
           }
           if(pth->energy_repart_coefficient==no_factorization) hyrec_data.cosmo->inj_params->energy_repart_coefficient = 0;
           else if(pth->energy_repart_coefficient==GSVI || pth->energy_repart_coefficient ==chi_from_file) hyrec_data.cosmo->inj_params->energy_repart_coefficient = 1;
           else if(pth->energy_repart_coefficient==SSCK) hyrec_data.cosmo->inj_params->energy_repart_coefficient = 2;
           hyrec_data.cosmo->inj_params->f_esc = pth->f_esc;
           hyrec_data.cosmo->inj_params->Zeta_ion = pth->Zeta_ion ; /**< Lyman continuum photon production efficiency of the stellar population */
           hyrec_data.cosmo->inj_params->fx = pth->fx; /**< X-ray efficiency fudge factor of photons responsible for heating the medium. */
           hyrec_data.cosmo->inj_params->Ex = pth->Ex*_eV_over_joules_; /**< Associated normalization from Pober et al. 1503.00045. */
           hyrec_data.cosmo->inj_params->ap = pth->ap;   /**<  a few parameters entering the fit of the star formation rate (SFR), introduced in Madau & Dickinson, Ann.Rev.Astron.Astrophys. 52 (2014) 415-486, updated in Robertson & al. 1502.02024.*/
           hyrec_data.cosmo->inj_params->bp = pth->bp;
           hyrec_data.cosmo->inj_params->cp = pth->cp;
           hyrec_data.cosmo->inj_params->dp = pth->dp;
           hyrec_data.cosmo->inj_params->z_start_reio_stars = pth->z_start_reio_stars; /**< Controls the beginning of star reionisation, the SFR experiences is put to 0 above this value. */
           hyrec_data.cosmo->fsR = alpha_ratio;
           hyrec_data.cosmo->meR = me_ratio;
           hyrec_data.cosmo->inj_params->Mpbh = pth->PBH_accreting_mass;
           hyrec_data.cosmo->inj_params->fpbh = pth->PBH_fraction;
           hyrec_data.cosmo->inj_params->coll_ion = pth->coll_ion_pbh;
           hyrec_data.cosmo->inj_params->PBH_ADAF_delta = pth->PBH_ADAF_delta;
           hyrec_data.cosmo->inj_params->PBH_accretion_eigenvalue = pth->PBH_accretion_eigenvalue;
           if(pth->PBH_accretion_recipe == spherical_accretion) hyrec_data.cosmo->inj_params->PBH_accretion_recipe = 0;
           else if(pth->PBH_accretion_recipe == disk_accretion) hyrec_data.cosmo->inj_params->PBH_accretion_recipe = 1;
          //

   if (pth->thermodynamics_verbose > 0)
     printf(" -> calling HyRec version %s,\n",HYREC_VERSION);
   if (pth->thermodynamics_verbose > 0)
     printf("by Y. Ali-Haïmoud & C. Hirata\n");

   hyrec_compute_CLASS(&hyrec_data, FULL);


   /** - fill a few parameters in preco and pth */

   Nz=ppr->recfast_Nz0;

   preco->rt_size = Nz;
   preco->H0 = pba->H0 * _c_ / _Mpc_over_m_;
   /* preco->H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
   preco->YHe = pth->YHe;
   preco->Nnow = 3.*preco->H0*preco->H0*pba->Omega0_b*(1.-preco->YHe)/(8.*_PI_*_G_*_m_H_);
   /* energy injection parameters */
   preco->annihilation = pth->annihilation;
   preco->has_on_the_spot = pth->has_on_the_spot;
   preco->decay_fraction = pth->decay_fraction;
   preco->annihilation_f_halo = pth->annihilation_f_halo;
   preco->annihilation_z_halo = pth->annihilation_z_halo;
   preco->has_UCMH_spike = pth->has_UCMH_spike; // GFA
   preco->boost_table = pth->boost_table;
   preco->z_table_for_boost = pth->z_table_for_boost;
   pth->n_e=preco->Nnow;

   /** - allocate memory for thermodynamics interpolation tables (size known in advance) and fill it */

   class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);

   for(i = 0; i < Nz; i++) {

     /** - --> get redshift, corresponding results from hyrec, and background quantities */

     z = ppr->recfast_z_initial * (1. - (double)(i+1) / (double)Nz);

     xe = hyrec_xe(z, &hyrec_data);
     Tm = hyrec_Tm(z, &hyrec_data);

     class_call(background_tau_of_z(pba,
                                    z,
                                    &tau),
                pba->error_message,
                pth->error_message);

     class_call(background_at_tau(pba,
                                  tau,
                                  pba->short_info,
                                  pba->inter_normal,
                                  &last_index_back,
                                  pvecback),
                pba->error_message,
                pth->error_message);

     /*   class_call(thermodynamics_energy_injection(ppr,pba,preco,z,&energy_rate,pth->error_message),
          pth->error_message,
          pth->error_message);
     */

     /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
     Hz=pvecback[pba->index_bg_H] * _c_ / _Mpc_over_m_;

     /** - --> store the results in the table */

     /* results are obtained in order of decreasing z, and stored in order of growing z */

     /* redshift */
     *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z)=z;

     /* ionization fraction */
     *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe)=xe;

     /* Tb */
     *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb)=Tm;

     /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)
        with (1+z)dlnTb/dz= - [dlnTb/dlna] */
     *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
       = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + xe * (1.-pth->YHe)) * Tm *(1. - hyrec_dTmdlna(z, &hyrec_data) / Tm / 3.);

     /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
     *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
       = (1.+z) * (1.+z) * preco->Nnow * xe * _sigma_ * _Mpc_over_m_;

   }

   /* Cleanup */
   free(buffer);
   hyrec_free(&hyrec_data);
 #else

   class_stop(pth->error_message,
              "you compiled without including the HyRec code, and now wish to use it. Either set the input parameter 'recombination' to something else than 'HyRec', or recompile after setting in the Makefile the appropriate path HYREC=... ");

 #endif
   return _SUCCESS_;
 }


int fill_recombination_structure(struct precision * ppr,
                                 struct background * pba,
                                 struct thermo * pth,
                                 struct recombination * preco){

   double mu_H,Lalpha,Lalpha_He,DeltaB,DeltaB_He;

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) */
  preco->rt_size = ppr->recfast_Nz0;
  class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);


  /** - read a few precision/cosmological parameters */


  /* preco->H0 is H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
  preco->H0 = pba->H0 * _c_ / _Mpc_over_m_;

  /* Yp */
  preco->YHe = pth->YHe;

  /* Tnow */
  preco->Tnow = pba->T_cmb;

  /* H_frac */
  preco->H_frac = ppr->recfast_H_frac;

  /* H fudging */
 class_test((ppr->recfast_Hswitch != _TRUE_) && (ppr->recfast_Hswitch != _FALSE_),
            pth->error_message,
            "RECFAST error: unknown H fudging scheme");
  preco->fu = ppr->recfast_fudge_H;
  if (ppr->recfast_Hswitch == _TRUE_)
    preco->fu += ppr->recfast_delta_fudge_H;

  /* He fudging */
  class_test((ppr->recfast_Heswitch < 0) || (ppr->recfast_Heswitch > 6),
             pth->error_message,
             "RECFAST error: unknown He fudging scheme");

  mu_H = 1./(1.-preco->YHe);
  Lalpha = 1./_L_H_alpha_;
  Lalpha_He = 1./_L_He_2p_;
  DeltaB = _h_P_*_c_*(_L_H_ion_-_L_H_alpha_);
  DeltaB_He = _h_P_*_c_*(_L_He1_ion_-_L_He_2s_);

  preco->fHe = preco->YHe/(_not4_ *(1.-preco->YHe)); /* recfast 1.4 */
  preco->Nnow = 3.*preco->H0*preco->H0*pba->Omega0_b/(8.*_PI_*_G_*mu_H*_m_H_);

  if (pth->has_extended_PBH_MassFunc == _TRUE_) { // GFA, allocate energy injections per each mass
   class_alloc(preco->energy_rate_at_mass,pth->num_PBH_accreting_mass*sizeof(double),pth->error_message);
  }

  /* energy injection parameters */
  preco->annihilation = pth->annihilation;
  preco->has_on_the_spot = pth->has_on_the_spot;
  preco->decay_fraction = pth->decay_fraction;
  preco->PBH_accreting_mass = pth->PBH_accreting_mass;
  preco->table_PBH_accreting_mass = pth->table_PBH_accreting_mass; // GFA
  preco->num_PBH_accreting_mass = pth->num_PBH_accreting_mass; // GFA
  preco->PBH_ADAF_delta = pth->PBH_ADAF_delta;
  preco->PBH_accretion_eigenvalue = pth->PBH_accretion_eigenvalue;
  preco->PBH_relative_velocities = pth->PBH_relative_velocities;
  preco->PBH_accretion_recipe = pth->PBH_accretion_recipe;
  preco->energy_deposition_function = pth->energy_deposition_function;
  preco->PBH_evaporating_mass = pth->PBH_evaporating_mass;
  preco->PBH_fraction = pth->PBH_fraction;

  preco->PBH_table_is_initialized = pth->PBH_table_is_initialized;
  preco->PBH_table_z = pth->PBH_table_z;
  preco->PBH_table_mass = pth->PBH_table_mass;
  preco->PBH_table_mass_dd = pth->PBH_table_mass_dd;
  preco->PBH_table_F = pth->PBH_table_F;
  preco->PBH_table_F_dd = pth->PBH_table_F_dd;

  preco->energy_repart_coefficient = pth->energy_repart_coefficient;
  preco->annihilation_f_halo = pth->annihilation_f_halo;
  preco->annihilation_z_halo = pth->annihilation_z_halo;
  preco->f_eff = pth->f_eff;
  preco->has_UCMH_spike = pth->has_UCMH_spike; //GFA
  preco->boost_table = pth->boost_table;
  preco->z_table_for_boost = pth->z_table_for_boost;


  /* quantities related to constants defined in thermodynamics.h */
  //n = preco->Nnow * pow((1.+z),3);

  preco->CDB = DeltaB/_k_B_;
  preco->CDB_He = DeltaB_He/_k_B_;
  preco->CB1 = _h_P_*_c_*_L_H_ion_/_k_B_;
  preco->CB1_He1 = _h_P_*_c_*_L_He1_ion_/_k_B_;
  preco->CB1_He2 = _h_P_*_c_*_L_He2_ion_/_k_B_;
  preco->CR = 2.*_PI_*(_m_e_/_h_P_)*(_k_B_/_h_P_);
  preco->CK = pow(Lalpha,3)/(8.*_PI_);
  preco->CK_He = pow(Lalpha_He,3)/(8.*_PI_);
  preco->CL = _c_*_h_P_/(_k_B_*Lalpha);
  preco->CL_He = _c_*_h_P_/(_k_B_/_L_He_2s_);
  preco->CT = (8./3.) * (_sigma_/(_m_e_*_c_)) *
    (8.*pow(_PI_,5)*pow(_k_B_,4)/ 15./ pow(_h_P_,3)/pow(_c_,3));

  preco->Bfact = _h_P_*_c_*(_L_He_2p_-_L_He_2s_)/_k_B_;

  return _SUCCESS_;
}

/**
 * Integrate thermodynamics with RECFAST.
 *
 * Integrate thermodynamics with RECFAST, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_recombination, from thermodynamics_init().
 *
 *
 *******************************************************************************
 * RECFAST is an integrator for Cosmic Recombination of Hydrogen and Helium,
 * developed by Douglas Scott (dscott@astro.ubc.ca)
 * based on calculations in the paper Seager, Sasselov & Scott
 * (ApJ, 523, L1, 1999).
 * and "fudge" updates in Wong, Moss & Scott (2008).
 *
 * Permission to use, copy, modify and distribute without fee or royalty at
 * any tier, this software and its documentation, for any purpose and without
 * fee or royalty is hereby granted, provided that you agree to comply with
 * the following copyright notice and statements, including the disclaimer,
 * and that the same appear on ALL copies of the software and documentation,
 * including modifications that you make for internal use or for distribution:
 *
 * Copyright 1999-2010 by University of British Columbia.  All rights reserved.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO
 * REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
 * BY WAY OF EXAMPLE, BUT NOT LIMITATION,
 * U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF
 * MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
 * THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
 * ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
 *******************************************************************************
 *
 * Version 1.5: includes extra fitting function from
 *              Rubino-Martin et al. arXiv:0910.4383v1 [astro-ph.CO]
 *
 * @param ppr      Input: pointer to precision structure
 * @param pba      Input: pointer to background structure
 * @param pth      Input: pointer to thermodynamics structure
 * @param preco    Output: pointer to recombination structure
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 * @return the error status
 */

int thermodynamics_recombination_with_recfast(
                                              struct precision * ppr,
                                              struct background * pba,
                                              struct thermo * pth,
                                              struct recombination * preco,
                                              double * pvecback
                                              ) {

  /** Summary: */

  /** - define local variables */

  /* vector of variables to be integrated: x_H, x_He, Tmat */
  double y[3],dy[3];

  /* other recfast variables */
  double OmegaB,zinitial,x_He0,x0;
  double x_H0=0.;
  double z,mu_H,Lalpha,Lalpha_He,DeltaB,DeltaB_He;
  double zstart,zend,rhs;
  int i,Nz;

  /* introduced by JL for smoothing the various steps */
  double x0_previous,x0_new,s,weight;

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;

  /* contains all fixed parameters which should be passed to thermodynamics_derivs_with_recfast */
  struct thermodynamics_parameters_and_workspace tpaw;

  /** - Fill the recombination structure with all important parameters */
  class_call(fill_recombination_structure(ppr,pba,pth,preco),pth->error_message,pth->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */
  class_call(initialize_generic_integrator(_RECFAST_INTEG_SIZE_, &gi),
             gi.error_message,
             pth->error_message);

  /* Nz */
  Nz=ppr->recfast_Nz0;

  /* Omega_b */
  OmegaB = pba->Omega0_b;

  /* z_initial */
  zinitial=ppr->recfast_z_initial;

  /* related quantities */
  z=zinitial;
  // mu_H = 1./(1.-preco->YHe);
  pth->n_e = preco->Nnow;

  /* quantities related to constants defined in thermodynamics.h */

  // Lalpha = 1./_L_H_alpha_;
  // Lalpha_He = 1./_L_He_2p_;
  // DeltaB = _h_P_*_c_*(_L_H_ion_-_L_H_alpha_);
  // DeltaB_He = _h_P_*_c_*(_L_He1_ion_-_L_He_2s_);

  /** - define the fields of the 'thermodynamics parameter and workspace' structure */
  tpaw.pba = pba;
  tpaw.ppr = ppr;
  tpaw.preco = preco;
  tpaw.pth = pth;
  tpaw.pvecback = pvecback;

  /** - impose initial conditions at early times */

  class_test(zinitial < ppr->recfast_z_He_3,
             pth->error_message,
             "increase zinitial, otherwise should get initial conditions from recfast's get_init routine (less precise anyway)");

  y[0] = 1.;
  y[1] = 1.;
  x0 = 1.+2.*preco->fHe;
  y[2] = preco->Tnow*(1.+z);

  /** - loop over redshift steps Nz; integrate over each step with
      generic_integrator(), store the results in the table using
      thermodynamics_derivs_with_recfast()*/

  for(i=0; i <Nz; i++) {

    zstart = zinitial * (double)(Nz-i) / (double)Nz;
    zend   = zinitial * (double)(Nz-i-1) / (double)Nz;

    z = zend;

    /** - --> first approximation: H and Helium fully ionized */

    if (z > ppr->recfast_z_He_1+ppr->recfast_delta_z_He_1) {
      x_H0 = 1.;
      x_He0 = 1.;
      x0 = 1.+2.*preco->fHe;
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** - --> second approximation: first Helium recombination (analytic approximation) */

    else if (z > ppr->recfast_z_He_2+ppr->recfast_delta_z_He_2) {
      x_H0 = 1.;
      x_He0 = 1.;

      rhs = exp( 1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He2/(preco->Tnow*(1.+z)) ) / preco->Nnow;

      /* smoothed transition */
      if (z > ppr->recfast_z_He_1-ppr->recfast_delta_z_He_1) {
        x0_previous = 1.+2.*preco->fHe;
        x0_new = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));

        /* get s from -1 to 1 */
        s = (ppr->recfast_z_He_1-z)/ppr->recfast_delta_z_He_1;
        /* infer f1(s) = smooth function interpolating from 0 to 1 */
        weight = f1(s);

        x0 = weight*x0_new+(1.-weight)*x0_previous;
      }
      /* transition finished */
      else {
        x0 = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));
      }

      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** - --> third approximation: first Helium recombination completed */

    else if (z > ppr->recfast_z_He_3+ppr->recfast_delta_z_He_3) {
      x_H0 = 1.;
      x_He0 = 1.;

      /* smoothed transition */
      if (z > ppr->recfast_z_He_2-ppr->recfast_delta_z_He_2) {
        rhs = exp( 1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He2/(preco->Tnow*(1.+z)) ) / preco->Nnow;
        x0_previous = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));
        x0_new = 1. + preco->fHe;
        /* get s from -1 to 1 */
        s = (ppr->recfast_z_He_2-z)/ppr->recfast_delta_z_He_2;
        /* infer f1(s) = smooth function interpolating from 0 to 1 */
        weight = f1(s);

        x0 = weight*x0_new+(1.-weight)*x0_previous;

      }
      /* transition finished */
      else {
        x0 = 1.+preco->fHe;
      }

      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** - --> fourth approximation: second Helium recombination starts (analytic approximation) */

    else if (y[1] > ppr->recfast_x_He0_trigger ) {
      x_H0 = 1.;

      rhs = 4.*exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He1/(preco->Tnow*(1.+z)))/preco->Nnow;
      x_He0 = 0.5*(sqrt(pow((rhs-1.),2) + 4.*(1.+preco->fHe)*rhs )- (rhs-1.));

      /* smoothed transition */
      if (z > ppr->recfast_z_He_3-ppr->recfast_delta_z_He_3) {
        x0_previous = 1. + preco->fHe;
        x0_new = x_He0;
        /* get s from -1 to 1 */
        s = (ppr->recfast_z_He_3-z)/ppr->recfast_delta_z_He_3;
        /* infer f1(x) = smooth function interpolating from 0 to 1 */
        weight = f1(s);

        x0 = weight*x0_new+(1.-weight)*x0_previous;
      }
      /* transition finished */
      else {
        x0 = x_He0;
      }

      x_He0 = (x0-1.)/preco->fHe;
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** - --> fifth approximation: second Helium recombination (full
        evolution for Helium), H recombination starts (analytic
        approximation) */

    else if (y[0] > ppr->recfast_x_H0_trigger && z > 200) {

      rhs = exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1/(preco->Tnow*(1.+z)))/preco->Nnow;
      x_H0 = 0.5*(sqrt(pow(rhs,2)+4.*rhs) - rhs);

      class_call(generic_integrator(thermodynamics_derivs_with_recfast,
                                    zstart,
                                    zend,
                                    y,
                                    &tpaw,
                                    ppr->tol_thermo_integration,
                                    ppr->smallest_allowed_variation,
                                    &gi),
                 gi.error_message,
                 pth->error_message);

      y[0] = MIN(x_H0,1);   //security added by VP: y[0] is only the hydrogen contribution and cannot be greater than 1

      if (pth->thermodynamics_verbose > 2) {
        fprintf(stdout, "in function thermodynamics_recombination_with_recfast, fifth approximation : zend %e y[0] %e\n",zend, y[0]);
      }
      /* smoothed transition */
      if (ppr->recfast_x_He0_trigger - y[1] < ppr->recfast_x_He0_trigger_delta && z > 200) {
        rhs = 4.*exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He1/(preco->Tnow*(1.+z)))/preco->Nnow;
        x0_previous = 0.5*(sqrt(pow((rhs-1.),2) + 4.*(1.+preco->fHe)*rhs )- (rhs-1.));
        x0_new = y[0] + preco->fHe*y[1];
        /* get s from 0 to 1 */
        s = (ppr->recfast_x_He0_trigger - y[1])/ppr->recfast_x_He0_trigger_delta;
        /* infer f2(x) = smooth function interpolating from 0 to 1 */
        weight = f2(s);

        x0 = weight*x0_new+(1.-weight)*x0_previous;
      }
      /* transition finished */
      else {
        x0 = y[0] + preco->fHe*y[1];
      }
      // x0 = y[0] + preco->fHe*y[1];

    }

    /** - --> last case: full evolution for H and Helium */

    else {

      /* quantities used for smoothed transition */
      if (ppr->recfast_x_H0_trigger - y[0] < ppr->recfast_x_H0_trigger_delta  && z > 200 ) {
        rhs = exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1/(preco->Tnow*(1.+z)))/preco->Nnow;
        x_H0 = 0.5*(sqrt(pow(rhs,2)+4.*rhs) - rhs);
      }
      else x_H0 = y[0];

      class_call(generic_integrator(thermodynamics_derivs_with_recfast,
                                    zstart,
                                    zend,
                                    y,
                                    &tpaw,
                                    ppr->tol_thermo_integration,
                                    ppr->smallest_allowed_variation,
                                    &gi),
                 gi.error_message,
                 pth->error_message);
       y[0] = MIN(y[0],1); //security added by VP: y[0] is only the hydrogen contribution and cannot be greater than 1


      /* smoothed transition */
      if (ppr->recfast_x_H0_trigger - y[0] < ppr->recfast_x_H0_trigger_delta && z > 200) {
        /* get s from 0 to 1 */
        s = (ppr->recfast_x_H0_trigger - y[0])/ppr->recfast_x_H0_trigger_delta;
        /* infer f2(s) = smooth function interpolating from 0 to 1 */
        weight = f2(s);

        x0 = weight*y[0]+(1.-weight)*x_H0 + preco->fHe*y[1];

      }
      /* transition finished */
      else {
        x0 = y[0] + preco->fHe*y[1];
      }

        // x0 = y[0] + preco->fHe*y[1];
        if(pth->thermodynamics_verbose>2){
          fprintf(stdout, "in function thermodynamics_recombination_with_recfast, full calculation zend %e x0 %e y[0] %e\n",zend, x0, y[0]);
        }

    }
  /*  double argument;
    argument = (pth->helium_fullreio_redshift - z)
        /pth->helium_fullreio_width;
      x0 += pth->YHe/(_not4_*(1.-pth->YHe))* (tanh(argument)+1.);*/

    /** - --> store the results in the table */
    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* redshift */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z)=zend;

    /* ionization fraction */

    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe)=x0;

    /* Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb)=y[2];

    /* get dTb/dz=dy[2] */
    class_call(thermodynamics_derivs_with_recfast(zend, y, dy, &tpaw,pth->error_message),
               pth->error_message,
               pth->error_message);

    /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
      = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * preco->YHe + x0 * (1.-preco->YHe)) * y[2] * (1. + (1.+zend) * dy[2] / y[2] / 3.);

    /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
      = (1.+zend) * (1.+zend) * preco->Nnow * x0 * _sigma_ * _Mpc_over_m_;
      if(pth->thermodynamics_verbose>1){
        fprintf(stdout,"%e %e %e %e %e %e\n",
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z),
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe),
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb),
             (1.+zend) * dy[2],
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2),
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
             );
      }


  }

  /** - cleanup generic integrator with cleanup_generic_integrator() */

  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             pth->error_message);

  return _SUCCESS_;
}

/**
 * Subroutine evaluating the derivative with respect to redshift of
 * thermodynamical quantities (from RECFAST version 1.4).
 *
 * Computes derivatives of the three variables to integrate: \f$ d x_H
 * / dz, d x_{He} / dz, d T_{mat} / dz \f$.
 *
 * This is one of the few functions in the code which are passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed parameters and workspaces are passed through a generic
 *   pointer. Here, this pointer contains the precision, background
 *   and recombination structures, plus a background vector, but
 *   generic_integrator() doesn't know its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 *   usual to pth->error_message, but to a generic error_message
 *   passed in the list of arguments.
 *
 * @param z                        Input: redshift
 * @param y                        Input: vector of variable to integrate
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices) and workspace (already allocated)
 * @param error_message            Output: error message
 */

int thermodynamics_derivs_with_recfast(
                                       double z,
                                       double * y,
                                       double * dy,
                                       void * parameters_and_workspace,
                                       ErrorMsg error_message
                                       ) {


  /* define local variables */
  double x,n,n_He,Trad,Tmat,x_H,x_He,Hz,dHdz,epsilon;
  double Rup,Rup_2,Rdown,K,K_He,Rup_He,Rup_He_2,Rdown_He,He_Boltz;
  double timeTh,timeH;
  double sq_0,sq_1;

  /* new in recfast 1.4: */
  double Rdown_trip,Rup_trip,tauHe_s,pHe_s,Doppler,gamma_2Ps,pb,qb,AHcon;
  double tauHe_t,pHe_t,CL_PSt,gamma_2Pt;
  double CfHe_t=0.;
  int Heflag;

  struct thermodynamics_parameters_and_workspace * ptpaw;
  struct precision * ppr;
  struct background * pba;
  struct thermo * pth;
  struct recombination * preco;
  double * pvecback;
  int last_index_back;

  /* used for energy injection from dark matter */
  double C;
  //double C_He;
  double energy_rate;
  // GFA
  double energy_rate_dep_ion;
  double energy_rate_dep_lya;
  double energy_rate_dep_heat;

  double tau;
  double chi_heat;
  double chi_lya;
  double chi_ionH;
  double chi_ionHe;
  double chi_lowE;
  double dTdz_DM, dTdz_CMB, dTdz_adia, dTdz_stars;
   /*used for reionization from realistic star model*/
  double rho_sfr,stars_xe,dNion_over_dt,L_x;

  int index_M; // GFA
  int i;
  double rho_cdm_today;
  double M_sun = 2e30; // in Kg

  ptpaw = parameters_and_workspace;
  ppr = ptpaw->ppr;
  pba = ptpaw->pba;
  pth = ptpaw->pth;
  preco = ptpaw->preco;
  pvecback = ptpaw->pvecback;
  rho_sfr = pth->ap*pow(1+z,pth->bp)/(1+pow((1+z)/pth->cp,pth->dp))/pow(_Mpc_over_m_,3)*pow(1+z,3)*(1+tanh((pth->z_start_reio_stars-z)))/2;//add a (sharp) smoothing function.

  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */

  /* security added by Vivian Poulin to avoid bug when dealing with energy injection modifying ionisation history */
  x_H = MIN(y[0],1.);
  x_H = MAX(y[0],0.);
  x_He = MIN(y[1],1.);
  x_He = MAX(y[1],0.);
  x = MIN(x_H + preco->fHe * x_He,1+preco->fHe);
  x = MAX(x_H + preco->fHe * x_He,0.);

  // x_H = y[0];
  // x_He = y[1];
  // x = x_H + preco->fHe * x_He;

  Tmat = MAX(y[2],0.);

  // fprintf(stderr, "input xH %e xHe %e x %e Tmat %e\n",x_H, x_He, x, Tmat );
  n = preco->Nnow * (1.+z) * (1.+z) * (1.+z);
  n_He = preco->fHe * n;
  Trad = preco->Tnow * (1.+z);
  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->short_info,
                               pba->inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             error_message);

   if((pth->annihilation!=0 || pth->decay_fraction!=0 || pth->PBH_accreting_mass!=0 || pth->PBH_evaporating_mass != 0 || pth->has_extended_PBH_MassFunc == _TRUE_)){ // GFA
     preco->xe_tmp=x;
     preco->Tm_tmp=Tmat;

        if( z > 2){//sometimes problem with interpolation

          if (pth->has_extended_PBH_MassFunc == _TRUE_) { // GFA, compute energy injection for each mass
            class_call(thermodynamics_accreting_pbh_energy_injection_PBH_MF(ppr,pba,preco,z,error_message),
                       error_message,
                       error_message);

          } else {
            class_call(thermodynamics_energy_injection(ppr,pba,preco,z,&energy_rate,error_message),
                       error_message,
                       error_message);
          }

        }
        else {
          energy_rate = 0;
          if (pth->has_extended_PBH_MassFunc == _TRUE_) {
          for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) {
          preco->energy_rate_at_mass[index_M] =  0.;
          }
          }

        }
    preco->z_tmp=z;
    } else {
      energy_rate=0;
    }

 // fprintf(stdout,"%e      %e     %e      %e      %e    \n", x,pth->chi_heat,pth->chi_lya, pth->chi_ionH,pth->chi_ionHe,pth->chi_lowE);
  /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
  Hz=pvecback[pba->index_bg_H]* _c_ / _Mpc_over_m_;
  //fprintf(stdout,"T_mat = %e\t_a_PPB_ = %e\t_b_PPB_ = %e\t_c_PPB_ = %e\t_d_PPB_ = %e\n",(double) Tmat,(double) _a_PPB_, (double) _b_PPB_, (double) _c_PPB_, (double) _d_PPB_);
  Rdown=1.e-19*_a_PPB_*pow((Tmat/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Tmat/1.e4),_d_PPB_));
  Rup_2 = 1.e-19*_a_PPB_*pow((Trad/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Trad/1.e4),_d_PPB_)) * pow((preco->CR*Trad),1.5)*exp(-preco->CDB/Trad);
  Rup = Rdown * pow((preco->CR*Tmat),1.5)*exp(-preco->CDB/Tmat);

  sq_0 = sqrt(Tmat/_T_0_);
  sq_1 = sqrt(Tmat/_T_1_);
  Rdown_He = _a_VF_/(sq_0 * pow((1.+sq_0),(1.-_b_VF_)) * pow((1. + sq_1),(1. + _b_VF_)));
  Rup_He_2 = 4.*Rdown_He*pow((preco->CR*Trad),1.5)*exp(-preco->CDB_He/Trad);
  Rup_He = 4.*Rdown_He*pow((preco->CR*Tmat),1.5)*exp(-preco->CDB_He/Tmat);
  K = preco->CK/Hz;

  /* following is from recfast 1.5 */

  if (ppr->recfast_Hswitch == _TRUE_ )
    K *= 1.
      + ppr->recfast_AGauss1*exp(-pow((log(1.+z)-ppr->recfast_zGauss1)/ppr->recfast_wGauss1,2))
      + ppr->recfast_AGauss2*exp(-pow((log(1.+z)-ppr->recfast_zGauss2)/ppr->recfast_wGauss2,2));

  /* end of new recfast 1.5 piece */

  /* following is from recfast 1.4 */

  Rdown_trip = _a_trip_/(sq_0*pow((1.+sq_0),(1.-_b_trip_)) * pow((1.+sq_1),(1.+_b_trip_)));
  Rup_trip = Rdown_trip*exp(-_h_P_*_c_*_L_He2St_ion_/(_k_B_*Tmat))*pow(preco->CR*Tmat,1.5)*4./3.;

  if ((x_He < 5.e-9) || (x_He > ppr->recfast_x_He0_trigger2))
    Heflag = 0;
  else
    Heflag = ppr->recfast_Heswitch;

  if (Heflag == 0)
    K_He = preco->CK_He/Hz;
  else {
    tauHe_s = _A2P_s_*preco->CK_He*3.*n_He*(1.-x_He)/Hz;
    pHe_s = (1.-exp(-tauHe_s))/tauHe_s;
    K_He = 1./(_A2P_s_*pHe_s*3.*n_He*(1.-x_He));

    /*    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.99999)) { */
    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.9999999)) { /* threshold changed by Antony Lewis in 2008 to get smoother Helium */

      Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
      Doppler = _c_*_L_He_2p_*sqrt(Doppler);
      gamma_2Ps = 3.*_A2P_s_*preco->fHe*(1.-x_He)*_c_*_c_
        /(sqrt(_PI_)*_sigma_He_2Ps_*8.*_PI_*Doppler*(1.-x_H))
        /pow(_c_*_L_He_2p_,2);
      pb = 0.36;
      qb = ppr->recfast_fudge_He;
      AHcon = _A2P_s_/(1.+pb*pow(gamma_2Ps,qb));
      K_He=1./((_A2P_s_*pHe_s+AHcon)*3.*n_He*(1.-x_He));
    }

    if (Heflag >= 3) {
      tauHe_t = _A2P_t_*n_He*(1.-x_He)*3./(8.*_PI_*Hz*pow(_L_He_2Pt_,3));
      pHe_t = (1. - exp(-tauHe_t))/tauHe_t;
      CL_PSt = _h_P_*_c_*(_L_He_2Pt_ - _L_He_2St_)/_k_B_;
      if ((Heflag == 3) || (Heflag == 5) || (x_H >= 0.99999)) {
        CfHe_t = _A2P_t_*pHe_t*exp(-CL_PSt/Tmat);
        CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
      else {
        Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
        Doppler = _c_*_L_He_2Pt_*sqrt(Doppler);
        gamma_2Pt = 3.*_A2P_t_*preco->fHe*(1.-x_He)*_c_*_c_
          /(sqrt(_PI_)*_sigma_He_2Pt_*8.*_PI_*Doppler*(1.-x_H))
          /pow(_c_*_L_He_2Pt_,2);
        pb = 0.66;
        qb = 0.9;
        AHcon = _A2P_t_/(1.+pb*pow(gamma_2Pt,qb))/3.;
        CfHe_t = (_A2P_t_*pHe_t+AHcon)*exp(-CL_PSt/Tmat);
        CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
    }
  }

  /* end of new recfast 1.4 piece */

  timeTh=(1./(preco->CT*pow(Trad,4)))*(1.+x+preco->fHe)/x;
  timeH=2./(3.*preco->H0*pow(1.+z,1.5));

  /************/
  /* hydrogen */
  /************/

  if (x_H > ppr->recfast_x_H0_trigger)
    dy[0] = 0.;
  else {
    /* equations modified to take into account energy injection from dark matter */
      chi_ionH = 0.;
      chi_ionHe = 0.;
      chi_lya = 0.;

    if(preco->annihilation > 0 || preco->decay_fraction > 0 || preco->PBH_accreting_mass > 0 || preco->PBH_evaporating_mass > 0 || pth->has_extended_PBH_MassFunc == _TRUE_){
      if (x < 1.){
        /* coefficient as revised by Galli et al. 2013 (in fact it is an interpolation by Vivian Poulin of Table V of Galli et al. 2013) */
        if(pth->energy_repart_coefficient==GSVI|| pth->energy_repart_coefficient ==chi_from_file){
          class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,x),
                         error_message,
                         error_message);
          chi_ionH = pth->chi_ionH;
          chi_ionHe = pth->chi_ionHe;
          chi_lya = pth->chi_lya;

        }
        if(pth->energy_repart_coefficient==no_factorization){
          if (pth->has_extended_PBH_MassFunc == _TRUE_){ // GFA, interpolate functions f(z) per channel and now also per each mass
            class_call(thermodynamics_annihilation_coefficients_interpolate_PBH_MF(ppr,pba,pth,z),
                           error_message,
                           error_message);
          } else {
            class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,z),
                           error_message,
                           error_message);
            chi_ionH = pth->chi_ionH;
            chi_ionHe = pth->chi_ionHe;
            chi_lya = pth->chi_lya;
          }

        }
        /* old approximation from Chen and Kamionkowski */
        if(pth->energy_repart_coefficient==SSCK){
          chi_ionH = (1.-x)/3.;
          chi_lya = chi_ionH;
          chi_ionHe=0;
        }

        chi_ionH = MIN(chi_ionH,1.);
        chi_ionHe = MIN(chi_ionHe,1.);
        chi_lya = MIN(chi_lya,1.);
        chi_ionH = MAX(chi_ionH,0.);
        chi_ionHe = MAX(chi_ionHe,0.);
        chi_lya = MAX(chi_lya,0.);
        if (pth->has_extended_PBH_MassFunc == _TRUE_) {
        for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) {
        pth->chi_ionH_at_mass[index_M] = MIN(pth->chi_ionH_at_mass[index_M],1.);
        pth->chi_ionHe_at_mass[index_M] = MIN(pth->chi_ionHe_at_mass[index_M],1.);
        pth->chi_lya_at_mass[index_M] = MIN(pth->chi_lya_at_mass[index_M],1.);
        pth->chi_ionH_at_mass[index_M] = MAX(pth->chi_ionH_at_mass[index_M],0.);
        pth->chi_ionHe_at_mass[index_M] = MAX(pth->chi_ionHe_at_mass[index_M],0.);
        pth->chi_lya_at_mass[index_M] = MAX(pth->chi_lya_at_mass[index_M],0.);
        }
        }

      }
      else {
        chi_ionH = 0.;
        chi_ionHe = 0.;
        chi_lya = 0.;

        if (pth->has_extended_PBH_MassFunc == _TRUE_) {
        for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) {
        pth->chi_ionH_at_mass[index_M] = 0.;
        pth->chi_ionHe_at_mass[index_M] = 0.;
        pth->chi_lya_at_mass[index_M] = 0.;
        }
        }

      }
      if(pth->thermodynamics_verbose>10){
        fprintf(stdout, "chi_ionH %e chi_ionHe %e chi_lya %e z% e\n", chi_ionH , chi_ionHe, chi_lya, z);
      }
    }

    /* Peebles' coefficient (approximated as one when the Hydrogen
           ionization fraction is very close to one) */
    if (x_H < ppr->recfast_x_H0_trigger2) {
      C = (1. + K*pth->Lambda_over_theoritical_Lambda*_Lambda_*n*(1.-x_H))/(1./preco->fu+K*pth->Lambda_over_theoritical_Lambda*_Lambda_*n*(1.-x_H)/preco->fu +K*Rup_2*n*(1.-x_H));  /* 2 modifications : 1) Rup -> Rup_2 evaluating the coefficient using Trad instead of Tmat; 2) add pth->Lambda_over_theoritical_Lambda, 1 in the standard case, allow to constraint A2s1s otherwise*/
      // C = (1. + K*_Lambda_*n*(1.-x_H))/(1./preco->fu+K*_Lambda_*n*(1.-x_H)/preco->fu +K*Rup*n*(1.-x_H));
      // fprintf(stdout, "A2s1s %e\n", pth->Lambda_over_theoritical_Lambda*_Lambda_);
    }
    else {
      C = 1.;

      }

      if (pth->has_extended_PBH_MassFunc == _TRUE_) {
       for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) { // evaluate the integrand per each mass
         pth->ener_rate_dep_ion_per_mass[index_M] = pth->table_PBH_MassFunc[index_M]*(pth->chi_ionH_at_mass[index_M]+pth->chi_ionHe_at_mass[index_M])*preco->energy_rate_at_mass[index_M];
         pth->ener_rate_dep_lya_per_mass[index_M] = pth->table_PBH_MassFunc[index_M]*pth->chi_lya_at_mass[index_M]*preco->energy_rate_at_mass[index_M];
       }

       energy_rate_dep_ion = 0.;
       energy_rate_dep_lya = 0.;
       for (i=0; i < pth->num_PBH_accreting_mass-1; i++) { //integrate over the mass using trapezoidal integral, improve to consider Simpson, although with at least 50 mass bins it shouldn't matter
         energy_rate_dep_ion += 0.5*(pth->table_PBH_accreting_mass[i+1]-pth->table_PBH_accreting_mass[i])*(pth->ener_rate_dep_ion_per_mass[i]+pth->ener_rate_dep_ion_per_mass[i+1]);
         energy_rate_dep_lya += 0.5*(pth->table_PBH_accreting_mass[i+1]-pth->table_PBH_accreting_mass[i])*(pth->ener_rate_dep_lya_per_mass[i]+pth->ener_rate_dep_lya_per_mass[i+1]);
       }
       energy_rate_dep_ion *= preco->PBH_fraction*(rho_cdm_today/(M_sun*_c_*_c_))*pow(1+z,3);
       energy_rate_dep_lya *= preco->PBH_fraction*(rho_cdm_today/(M_sun*_c_*_c_))*pow(1+z,3);

     } else { // standard case, just one mass or one energy injection
        energy_rate_dep_ion = energy_rate*(chi_ionH+chi_ionHe);
        energy_rate_dep_lya = energy_rate*chi_lya;
      }

      /* evolution of hydrogen ionisation fraction: */

      dy[0] = (x*x_H*n*Rdown - Rup_2*(1.-x_H)*exp(-preco->CL/Tmat)) * C / (Hz*(1.+z))       /* Peeble's equation with fudged factors */
            -(energy_rate_dep_ion/(_L_H_ion_*n)+energy_rate_dep_lya*(1.-C)/(_L_H_alpha_*n))/(_h_P_*_c_*Hz*(1.+z)); /* energy injection (neglect fraction going to helium) */

      // dy[0] = -5.89e-5*sqrt(Tmat/1e4)*exp(-1.58e5/Tmat)*(1-x_H)*x/ (Hz*(1.+z)) * 1e-6;     // Collisional ionisation, taken from 1503.04827, last factor is for conversion cm^3->m^3


      // fprintf(stdout, "z %e Tmat %e collision %e DM  %e standard %e\n",z, Tmat, -5.89e-5*sqrt(Tmat/1e4)*exp(-1.58e5/Tmat)*(1-x_H)*x/ (Hz*(1.+z)) * 1e-6,-energy_rate/n*((chi_ionH+chi_ionHe)/_L_H_ion_+chi_lya*(1.-C)/_L_H_alpha_)/(_h_P_*_c_*Hz*(1.+z)),(x*x_H*n*Rdown - Rup_2*(1.-x_H)*exp(-preco->CL/Tmat)) * C / (Hz*(1.+z)));
      if(pth->reio_parametrization == reio_stars_sfr_source_term){
        dNion_over_dt=pth->f_esc*pth->Zeta_ion*rho_sfr;
        stars_xe=dNion_over_dt/(Hz*(1.+z)*n);
        dy[0] -= stars_xe*(1-x)/3;
        // fprintf(stdout, " %e  %e  %e %e %e %e  %e\n",rho_sfr,stars_xe, dNion_over_dt,Hz,n,(1-x)/3,z );
      }
    // JL: test for debugginf reio_inter
    //fprintf(stdout,"%e  %e  %e  %e\n",z,Tmat,K*_Lambda_*n,K*Rup*n);

      if(pth->thermodynamics_verbose>10){
      fprintf(stdout, "z %e Tmat %e  DM  %e standard %e stars %e \n",z, Tmat,-energy_rate/n*((chi_ionH+chi_ionHe)/_L_H_ion_+chi_lya*(1.-C)/_L_H_alpha_)/(_h_P_*_c_*Hz*(1.+z)),(x*x_H*n*Rdown - Rup_2*(1.-x_H)*exp(-preco->CL/Tmat)) * C / (Hz*(1.+z)),stars_xe);
      }
  }

  /************/
  /* helium   */
  /************/

  if (x_He < 1.e-15)
    dy[1]=0.;
  else {

    if (preco->Bfact/Tmat < 680.)
      He_Boltz=exp(preco->Bfact/Tmat);
    else
      He_Boltz=exp(680.);

    /* equations modified to take into account energy injection from dark matter */
    //C_He=(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz)/(1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz);

    dy[1] = ((x*x_He*n*Rdown_He - Rup_He_2*(1.-x_He)*exp(-preco->CL_He/Tmat))
             *(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz))
      /(Hz*(1+z)* (1. + K_He*(_Lambda_He_+Rup_He_2)*n_He*(1.-x_He)*He_Boltz)); /* in case of energy injection due to DM, we neglect the contribution to helium ionization */

    // dy[1] += -2.02e-9*sqrt(Tmat/1e4)*exp(-2.85e5/Tmat)*(1-x_He)*x / (Hz*(1.+z)) * 1e-6;     // Collisional ionisation, taken from 1503.04827, last factor is for conversion cm^3->m^3
      //
      // /*******************Helium**********************/
      // dxedlna+=stars_xe*param->fHe*(1+tanh((6-z1)/0.5));
      // if(z1<6)dxedlna+=stars_xe*param->fHe*(1+tanh((3.5-z1)/0.5));
      /***********************************************/
    /* following is from recfast 1.4 */
    /* this correction is not self-consistent when there is energy injection  from dark matter, and leads to nan's  at small redshift (unimportant when reionization takes over before that redshift) */

    if (Heflag >= 3)
      dy[1] = dy[1] +
        (x*x_He*n*Rdown_trip
         - (1.-x_He)*3.*Rup_trip*exp(-_h_P_*_c_*_L_He_2St_/(_k_B_*Tmat)))
        *CfHe_t/(Hz*(1.+z));

    /* end of new recfast 1.4 piece */

  }

  if (timeTh < preco->H_frac*timeH) {
    /*   dy[2]=Tmat/(1.+z); */
    /* v 1.5: like in camb, add here a smoothing term as suggested by Adam Moss */
    dHdz=-pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_H]/pba->a_today* _c_ / _Mpc_over_m_;
    epsilon = Hz * (1.+x+preco->fHe) / (preco->CT*pow(Trad,3)*x);
    dy[2] = preco->Tnow + epsilon*((1.+preco->fHe)/(1.+preco->fHe+x))*((dy[0]+preco->fHe*dy[1])/x)
      - epsilon* dHdz/Hz + 3.*epsilon/(1.+z) ;
  }
  else {
    /* equations modified to take into account energy injection from dark matter */

    if(pth->annihilation >0 || pth->decay_fraction > 0 || pth->PBH_accreting_mass > 0 || pth->PBH_evaporating_mass > 0 || pth->has_extended_PBH_MassFunc == _TRUE_){
      if (x < 1.){
        /* coefficient as revised by Galli et al. 2013 (in fact it is an interpolation by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013) */
        if(pth->energy_repart_coefficient==GSVI|| pth->energy_repart_coefficient ==chi_from_file){
          class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,x),
                        error_message,
                        error_message);
            chi_heat = pth->chi_heat;
        }
        if(pth->energy_repart_coefficient==no_factorization){
          if (pth->has_extended_PBH_MassFunc == _TRUE_){  // GFA, interpolate functions f(z) per channel and now also per each mass
            class_call(thermodynamics_annihilation_coefficients_interpolate_PBH_MF(ppr,pba,pth,z),
                          error_message,
                          error_message);
          } else {
            class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,z),
                          error_message,
                          error_message);
              chi_heat = pth->chi_heat;
          }


        }
        /* old approximation from Chen and Kamionkowski */
        if(pth->energy_repart_coefficient==SSCK){
            chi_heat = (1.+2.*x)/3.;
        }

      } else {

        chi_heat = 1.;
        if (pth->has_extended_PBH_MassFunc == _TRUE_) {
        for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) {
        pth->chi_heat_at_mass[index_M] = 1.;
        }
        }

      }
        chi_heat= MIN(chi_heat,1.);
        chi_heat = MAX(chi_heat,0.);
        if (pth->has_extended_PBH_MassFunc == _TRUE_) {
        for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) {
        pth->chi_heat_at_mass[index_M]= MIN(pth->chi_heat_at_mass[index_M],1.);
        pth->chi_heat_at_mass[index_M] = MAX(pth->chi_heat_at_mass[index_M],0.);
        }
        }

    } else {
      chi_heat = 0.;
    }
    dTdz_adia=2.*Tmat/(1.+z);
    dTdz_CMB = preco->CT * pow(Trad,4) * x / (1.+x+preco->fHe) * (Tmat-Trad) / (Hz*(1.+z));
    if (pth->has_extended_PBH_MassFunc == _TRUE_) {
     for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) { // evaluate integrand at each mass
       pth->ener_rate_dep_heat_per_mass[index_M] = pth->table_PBH_MassFunc[index_M]*pth->chi_heat_at_mass[index_M]*preco->energy_rate_at_mass[index_M];
     }
     energy_rate_dep_heat = 0.;
     for (i=0; i < pth->num_PBH_accreting_mass-1; i++) { //GFA, integrate over mass using trapezoidal rule, improve to simpson, although with at least 50 mass bins it shouldn't matter
       energy_rate_dep_heat += 0.5*(pth->table_PBH_accreting_mass[i+1]-pth->table_PBH_accreting_mass[i])*(pth->ener_rate_dep_heat_per_mass[i]+pth->ener_rate_dep_heat_per_mass[i+1]);
     }
     energy_rate_dep_heat *= preco->PBH_fraction*(rho_cdm_today/(M_sun*_c_*_c_))*pow(1+z,3);
   } else { // standard case with single mass or single energy injection
     energy_rate_dep_heat = energy_rate*chi_heat;
     }
    dTdz_DM = -2./(3.*_k_B_)*energy_rate_dep_heat/n/(1.+preco->fHe+x)/(Hz*(1.+z));
    if(pth->star_heating_parametrization == heating_stars_sfr_source_term){
    L_x = 2*pth->Ex * pth->fx * rho_sfr/(3*_k_B_*n*Hz*(1.+z)*(1.+x+preco->fHe));
    dTdz_stars = -L_x*(1+2*x)/3.;
    }
    else dTdz_stars = 0;
    dy[2]= dTdz_CMB + dTdz_adia + dTdz_DM + dTdz_stars; /* dTdz_DM = energy injection */

    if(pth->thermodynamics_verbose>10){
      fprintf(stdout, "chi_heat %e z% e\n", chi_heat, z);
      fprintf(stdout, "z %e dT %e Tmat %e Trad %e dTdz_adia %e dTdz_CMB %e dTdz_DM %e dTdz_stars %e x %e\n", z, dy[2], Tmat, Trad,dTdz_adia, dTdz_CMB ,dTdz_DM,dTdz_stars,x);
      }

    // dy[2] += 5.89e-5*sqrt(Tmat/1e4)*exp(-1.58e5/Tmat)*(1-x_H)*x/ (Hz*(1.+z)) * 1e-6;
  }



  return _SUCCESS_;
}

/**
 * This routine merges the two tables 'recombination_table' and
 * 'reionization_table' inside the table 'thermodynamics_table', and
 * frees the temporary structures 'recombination' and 'reionization'.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pth   Input/Output: pointer to thermo structure
 * @param preco Input: pointer to filled recombination structure
 * @param preio Input: pointer to reionization structure
 * @return the error status
 */

int thermodynamics_merge_reco_and_reio(
                                       struct precision * ppr,
                                       struct thermo * pth,
                                       struct recombination * preco,
                                       struct reionization * preio
                                       ) {
  /** Summary: */

  /** - define local variables */

  int i,index_th,index_re;
  int index_M;

  /** - first, a little check that the two tables match each other and can be merged */

  if ((pth->reio_parametrization != reio_none)) {
    class_test(preco->recombination_table[preio->index_reco_when_reio_start*preco->re_size+preco->index_re_z] !=
               preio->reionization_table[(preio->rt_size -1)*preio->re_size+preio->index_re_z],
               pth->error_message,
               "mismatch which should never happen");
  }

  /** - find number of redshift in full table = number in reco + number in reio - overlap */

  pth->tt_size = ppr->recfast_Nz0 + preio->rt_size - preio->index_reco_when_reio_start - 1;


  /** - allocate arrays in thermo structure */

  class_alloc(pth->z_table,pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->thermodynamics_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->d2thermodynamics_dz2_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);


  /** - fill these arrays */

  for (i=0; i < preio->rt_size; i++) {
    pth->z_table[i]=
      preio->reionization_table[i*preio->re_size+preio->index_re_z];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_xe]=
      preio->reionization_table[i*preio->re_size+preio->index_re_xe];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa]=
      preio->reionization_table[i*preio->re_size+preio->index_re_dkappadtau];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_Tb]=
      preio->reionization_table[i*preio->re_size+preio->index_re_Tb];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_cb2]=
      preio->reionization_table[i*preio->re_size+preio->index_re_cb2];
  }
  for (i=0; i < ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 1; i++) {
    index_th=i+preio->rt_size;
    index_re=i+preio->index_reco_when_reio_start+1;
    pth->z_table[index_th]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_z];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_xe]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_xe];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_dkappa]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_dkappadtau];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_Tb]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_Tb];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_cb2]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_cb2];
  }

  /** - free the temporary structures */

  free(preco->recombination_table);
  if(pth->annihilation >0 || pth->decay_fraction > 0 || pth->PBH_accreting_mass > 0 || pth->PBH_evaporating_mass > 0){
    if(pth->has_on_the_spot == _FALSE_ && pth->energy_repart_coefficient!=no_factorization){
      thermodynamics_annihilation_f_eff_free(preco);
    }
    if(pth->energy_repart_coefficient == GSVI || pth->energy_repart_coefficient==no_factorization || pth->energy_repart_coefficient ==chi_from_file){
      thermodynamics_annihilation_coefficients_free(pth);
    }
  }
  if ((pth->reio_parametrization != reio_none))
    free(preio->reionization_table);

  if ((preco->PBH_table_is_initialized == _TRUE_) && pth->PBH_evaporating_mass > 0.) {
    // fprintf(stdout,"PBH tabels are free'd\n");
    free(preco->PBH_table_z);
    free(preco->PBH_table_mass);
    free(preco->PBH_table_mass_dd);
    free(preco->PBH_table_F);
    free(preco->PBH_table_F_dd);
    pth->PBH_table_is_initialized == _FALSE_;
  }

  if (pth->has_extended_PBH_MassFunc == _TRUE_) {
    // GFA: free all the matrices and vectors that have been allocated
    for (index_M=0; index_M < pth->num_PBH_accreting_mass; index_M++) {
      free(pth->annihil_coef_xe_at_mass[index_M]);
      free(pth->annihil_coef_heat_at_mass[index_M]);
      free(pth->annihil_coef_lya_at_mass[index_M]);
      free(pth->annihil_coef_ionH_at_mass[index_M]);
      free(pth->annihil_coef_ionHe_at_mass[index_M]);
      free(pth->annihil_coef_lowE_at_mass[index_M]);
      free(pth->annihil_coef_dd_heat_at_mass[index_M]);
      free(pth->annihil_coef_dd_lya_at_mass[index_M]);
      free(pth->annihil_coef_dd_ionH_at_mass[index_M]);
      free(pth->annihil_coef_dd_ionHe_at_mass[index_M]);
      free(pth->annihil_coef_dd_lowE_at_mass[index_M]);
    }

    free(pth->annihil_coef_xe_at_mass);
    free(pth->annihil_coef_heat_at_mass);
    free(pth->annihil_coef_lya_at_mass);
    free(pth->annihil_coef_ionH_at_mass);
    free(pth->annihil_coef_ionHe_at_mass);
    free(pth->annihil_coef_lowE_at_mass);
    free(pth->annihil_coef_dd_heat_at_mass);
    free(pth->annihil_coef_dd_lya_at_mass);
    free(pth->annihil_coef_dd_ionH_at_mass);
    free(pth->annihil_coef_dd_ionHe_at_mass);
    free(pth->annihil_coef_dd_lowE_at_mass);

    free(pth->table_PBH_accreting_mass);
    free(pth->table_PBH_MassFunc);
    free(pth->chi_heat_at_mass);
    free(pth->chi_lya_at_mass);
    free(pth->chi_ionH_at_mass);
    free(pth->chi_ionHe_at_mass);
    free(pth->chi_lowE_at_mass);
    free(pth->annihil_coef_num_lines_at_mass);
    free(preco->energy_rate_at_mass);

    free(pth->ener_rate_dep_ion_per_mass);
    free(pth->ener_rate_dep_lya_per_mass);
    free(pth->ener_rate_dep_heat_per_mass);

  }

  return _SUCCESS_;
}

/**
 * Subroutine for formatting thermodynamics output
 */

int thermodynamics_output_titles(struct background * pba,
                                 struct thermo *pth,
                                 char titles[_MAXTITLESTRINGLENGTH_]
                                 ){

  class_store_columntitle(titles,"z",_TRUE_);
  class_store_columntitle(titles,"conf. time [Mpc]",_TRUE_);
  class_store_columntitle(titles,"x_e",_TRUE_);
  class_store_columntitle(titles,"kappa' [Mpc^-1]",_TRUE_);
  //class_store_columntitle(titles,"kappa''",_TRUE_);
  //class_store_columntitle(titles,"kappa'''",_TRUE_);
  class_store_columntitle(titles,"exp(-kappa)",_TRUE_);
  class_store_columntitle(titles,"g [Mpc^-1]",_TRUE_);
  //class_store_columntitle(titles,"g'",_TRUE_);
  //class_store_columntitle(titles,"g''",_TRUE_);
  class_store_columntitle(titles,"Tb [K]",_TRUE_);
  class_store_columntitle(titles,"c_b^2",_TRUE_);
  class_store_columntitle(titles,"tau_d",_TRUE_);
  //class_store_columntitle(titles,"max. rate",_TRUE_);
  class_store_columntitle(titles,"r_d",pth->compute_damping_scale);
  class_store_columntitle(titles,"boost",pth->has_UCMH_spike); //GFA


  return _SUCCESS_;
}

int thermodynamics_output_data(struct background * pba,
                               struct thermo *pth,
                               int number_of_titles,
                               double *data
                               ){

  int index_z, storeidx;
  double *dataptr, *pvecthermo;
  double z,tau;
  double boost_at_z;

  //  pth->number_of_thermodynamics_titles = get_number_of_titles(pth->thermodynamics_titles);
  //pth->size_thermodynamics_data = pth->number_of_thermodynamics_titles*pth->tt_size;


  /* Store quantities: */
  for (index_z=0; index_z<pth->tt_size; index_z++){
    dataptr = data + index_z*number_of_titles;
    pvecthermo = pth->thermodynamics_table+index_z*pth->th_size;
    z = pth->z_table[index_z];
    storeidx=0;

    class_call(background_tau_of_z(
                                   pba,
                                   z,
                                   &tau
                                   ),
               pba->error_message,
               pth->error_message);

    class_store_double(dataptr,z,_TRUE_,storeidx);
    class_store_double(dataptr,tau,_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_xe],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_dkappa],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_ddkappa],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_dddkappa],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_exp_m_kappa],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_g],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_dg],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_ddg],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_Tb],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_cb2],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_tau_d],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_rate],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_r_d],pth->compute_damping_scale,storeidx);

    if (pth->has_UCMH_spike == _TRUE_) { //GFA
      if (z>1.e-3) {
        boost_at_z = array_interpolate_linear_simpler(pth->z_table_for_boost,pth->Number_z,pth->boost_table,z); //check ppr->Number_z
      } else {
        boost_at_z = pth->boost_table[0];
      }
    class_store_double(dataptr,boost_at_z,_TRUE_,storeidx);
    }

  }

  return _SUCCESS_;
}

int thermodynamics_tanh(double x,
                        double center,
                        double before,
                        double after,
                        double width,
                        double * result) {

  *result = before + (after-before)*(tanh((x-center)/width)+1.)/2.;

  return _SUCCESS_;
}

//GFA
int compute_boost_NFW_UCMH(
                           struct precision *  ppr,
                           struct background * pba,
                           struct thermo * pth
                          ) {

 struct halos_workspace params;
 double M_step, Mass_thres, M;
 double z_step, z;
 int index_z;
 int index_M;
 double rho_m_0, H100=1./3000.; /* #To express Hubble constant in units of 100 Mpc^{-1} (having set c=1) */
 double nu_minus,nu_plus,omega_of_z, N_frac_spike, c_spike, rho_c_over_rho_m;
 double S_standard=0., S_spike=0., S_tot=0.;
 double boost_low_mass, boost_high_mass, boost_spike;
 double gamma = 6.*pow(_PI_,2);
 double zF_min, zF_max;

 params.ppr = ppr;
 params.pba = pba;
 params.pth = pth;

 rho_m_0 = (pba->Omega0_cdm+pba->Omega0_b)*2.775e11*pow(pba->h,2);  /* present matter density in units of M_sol/Mpc^3 */
 Mass_thres = gamma*pow(pth->k_spike,-3)*rho_m_0;
 pth->Mass_min = gamma*pow(pth->k_fs,-3)*rho_m_0; //We always fix a minimal mass different from zero, even when we don't consider a free-streaming suppression in the transfer function
 if (pth->add_suppression_kfs_UCMH == _TRUE_) {
   pth->Mass_min *= 1.e-6; //if we add a suppression in the transfer function, we can integrate the mass function from arbitrarily small masses,
                            // since the free-streaming cut-off in the transfer will automatically kill the boost for M < M_min
 }
 if (pth->thermodynamics_verbose > 0) {
   printf("Computing boost factor for spiky primordial spectrum \n");
 }


 class_test((Mass_thres < pth->Mass_min),
            pth->error_message,
            "Mass_thres =%e Msun is smaller than Mass_min =%e Msun, this should never happen",Mass_thres, pth->Mass_min);

 S_spike =  (4./25.)*pth->A_spike*pow(pow(pth->k_spike,2.)*D_growth(0.,pba)*Transfer_f(pth->k_spike, pba,pth)/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h*pow(H100,2.)),2.);
 S_standard = integrate_simpson(ppr->k_min, pth->k_spike ,ppr->Number_k, LOG, integrand_for_S, &params);
 S_tot = S_spike + S_standard;

 z_step=(log10(ppr->z_max)-log10(ppr->z_min))/(ppr->Number_z-1);

 for (index_z=0; index_z < ppr->Number_z; index_z++) {
  z = ppr->z_min*pow(10.,index_z*z_step);
  pth->z_table_for_boost[index_z] = z;
  params.z = z;
  if (pth->UCMH_recipe == Delos) {
    // compute spike contribution (the only present in this case)
    zF_min = z;
    zF_max = pba->z_eq; //We restrict to mini-halos forming after matter-radiation equality, although next integral shouldn't depend a lot on this choice
    boost_spike = integrate_simpson(zF_min, zF_max,1000, LOG, integrand_boost_spike,&params);
    boost_spike *= pow(1.+z,-3.);
    if (boost_spike <0.) {
      boost_spike=0.;
    }
    pth->boost_table[index_z] = 1.+boost_spike;

  } else if (pth->UCMH_recipe == GG) {

    if (pth->consider_only_spike_UCMH == _FALSE_) {
      //Compute high-mass contribution
      boost_high_mass  = integrate_simpson(Mass_thres, ppr->Mass_max ,ppr->Number_M, LOG, integrand_boost_high_mass,&params);
      //Compute low-mass contribution
      boost_low_mass = integrate_simpson(pth->Mass_min, Mass_thres ,ppr->Number_M, LOG, integrand_boost_low_mass,&params);
      // Compute spike contribution
      if (pth->consider_zF_avg_UCMH == _TRUE_) { pth->M_at_Mthres = YES; }  //this is used to tell when to compute the average redshift at M = Mass_thres (i.e. at the spike)
      omega_of_z = _delta_crit_*D_growth(0.,pba)/D_growth(z,pba);
      nu_minus = omega_of_z/sqrt(S_tot);
      nu_plus  = omega_of_z/sqrt(S_standard);
      N_frac_spike = erf(nu_plus/sqrt(2.)) - erf(nu_minus/sqrt(2.));
      c_spike  = c_UCMH(Mass_thres, &params);
      boost_spike = N_frac_spike*one_halo_boost_UCMH(c_spike, Mass_thres, &params);
      if (pth->consider_zF_avg_UCMH == _TRUE_) {pth->M_at_Mthres = NO;}
      // add up three contribution
      pth->boost_table[index_z] = 1.+boost_spike+boost_low_mass+boost_high_mass;
      if (pth->UCMH_DM_ann_type == p_wave) {pth->boost_table[index_z] -= 1.;} //This is done because for p-wave annihilations we neglect the contribution from the smooth background (which is extremely small)

    } else {
      // Compute only spike contribution
      if (pth->consider_zF_avg_UCMH == _TRUE_) { pth->M_at_Mthres = YES;}  //this is used to tell when to compute the average redshift at M = Mass_thres (i.e. at the spike)
      omega_of_z = _delta_crit_*D_growth(0.,pba)/D_growth(z,pba);
      nu_minus = omega_of_z/sqrt(S_spike);
      N_frac_spike = erfc(nu_minus/sqrt(2.));
      c_spike  = c_UCMH(Mass_thres, &params);
      boost_spike = N_frac_spike*one_halo_boost_UCMH(c_spike, Mass_thres, &params);
      if (pth->consider_zF_avg_UCMH == _TRUE_) {pth->M_at_Mthres = NO;}
      pth->boost_table[index_z] = 1.+boost_spike;
      if (pth->UCMH_DM_ann_type == p_wave) {pth->boost_table[index_z] -= 1.;} //This is done because for p-wave annihilations we neglect the contribution from the smooth background (which is extremely small)
    }
  }

 }

 return _SUCCESS_;
}


/* required for CDM piece of transfer function, defined below */
double T0_tilde(double k,
                double alpha,
                double beta,
                struct background * pba) {
double C;
double q = k/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h); /* k needs to be passed in Mpc^{-1} */
C = (14.2/alpha) + 386./(1.+69.9*pow(q,1.08));
return log(exp(1.)+1.8*beta*q)/(log(exp(1.)+1.8*beta*q)+C*pow(q,2));
}

double T_Hu_no_baryon(double k,
                      struct background * pba) {
double C0, L0;
double q = k/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h); /* k needs to be passed in Mpc^{-1} */
L0 = log(2.*exp(1.)+1.8*q);
C0 = 14.2 + 731./(1.+62.5*q);
return L0/(L0+C0*pow(q,2));
}

/* required for baryon piece of transfer function, defined below */
double G(double y) {
return y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)));
}

/* define TRANSFER FUNCTION as given by Eisenstein-Hu 1998 (arXiv: 9709112) */
double T_Hu(double k,
            struct background * pba) {
/* k needs to be passed in Mpc^{-1} */
double z_eq, k_eq, z_d;
double b1, b2, b11, b22, a1, a2, alpha_c, beta_c;
double T_b, T_c, f;
double R_d, R_eq, s, k_silk;
double beta_node, alpha_b, beta_b, s_tilde;
double omega_m, omega_b;
omega_m = (pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h;
omega_b = pba->Omega0_b*pba->h*pba->h;
z_eq = (2.5e4)*omega_m;
k_eq = (7.46e-2)*omega_m;
b1 = 0.313*pow(omega_m,-0.419)*(1.+0.607*pow(omega_m,0.674));
b2 = 0.238*pow(omega_m,0.223);
z_d = 1291*(pow(omega_m,0.251)/(1.+0.659*pow(omega_m,0.828)))*(1.+b1*pow(omega_b,b2));
R_d = 31.5*omega_b*pow(z_d/1000.,-1);
R_eq = 31.5*omega_b*pow(z_eq/1000.,-1);
s = (2./(3.*k_eq))*sqrt(6./R_eq)*log((sqrt(1.+R_d)+sqrt(R_d+R_eq))/(1.+sqrt(R_eq)));
k_silk = 1.6*pow(omega_b,0.52)*pow(omega_m,0.73)*(1.+pow(10.4*omega_m,-0.95)); // in units of Mpc^{-1}
b11 = 0.944*pow(1.+ pow(458.*omega_m, -0.708),-1);
b22 = pow(0.395*omega_m,-0.0266);
beta_c = pow(1.+b11*(pow((omega_m-omega_b)/omega_m,b22)-1.),-1);
a1 = pow(46.9*omega_m, 0.670)*(1.+pow(32.1*omega_m, -0.532));
a2 = pow(12.0*omega_m, 0.424)*(1.+pow(45.0*omega_m, -0.582));
alpha_c  = pow(a1,-omega_b/omega_m)*pow(a2,-pow(omega_b/omega_m,3));
f = 1./(1.+pow(k*s/5.4,4));
T_c  = f*T0_tilde(k,1.,beta_c,pba) + (1.-f)*T0_tilde(k,alpha_c,beta_c,pba);
beta_node = 8.41*pow(omega_m,0.435);
s_tilde = s*pow(1.+pow(beta_node/(k*s),3),-1./3.);
beta_b = 0.5+(omega_b/omega_m)+(3.-2.*omega_b/omega_m)*sqrt(pow(17.2*omega_m,2)+1.);
alpha_b = 2.07*k_eq*s*pow(1.+R_d, -3./4.)*G((1.+z_eq)/(1.+z_d));
T_b = (T0_tilde(k,1.,1.,pba)/(1.+pow(k*s/5.2,2)) + (alpha_b/(1.+pow(beta_b/(k*s),3)))*exp(-pow(k/k_silk,1.4)))*sin(k*s_tilde)/(k*s_tilde);
return (omega_b/omega_m)*T_b + ((omega_m-omega_b)/omega_m)*T_c;
}


double integrand_for_S(void * params,
                            double k) { /* Notice we are not considering the z-dependence here, it has been factored-out */
 struct halos_workspace * params_local;
 struct background * pba;
 struct thermo * pth;
 double H100=1./3000.; /* #To express Hubble constant in units of 100 Mpc^{-1} (having set c=1) */
 double transfer, primordial;
 double cut_off = 1.;
 params_local = params;
 pba = params_local->pba;
 pth = params_local->pth;
 /* Ratio of matter spectrum and primordial spectrum, as given in Dodelson's book */
  transfer = (4.0/25.0)*pow(pow(k,2)*D_growth(0.,pba)*Transfer_f(k, pba,pth)/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h*pow(H100,2)),2);
 /* select the primordial spectrum  */
  primordial = pth->A_s*pow(cut_off,2)*pow(k/pth->k_pivot ,pth->n_s-1.);

 return primordial*transfer/k;

}


double integrand_boost_high_mass(void * params,
                                 double M) {
  struct halos_workspace * params_local;
  double c, integrand;
  params_local = params;
  c = c_NFW(M, params_local);
  integrand = M*halo_function_high_mass(params_local, M)*one_halo_boost_NFW(c, M, params_local);
  return integrand;
}

/* define function of concentration appearing for NFW profiles */
double one_halo_boost_NFW(double c,
                          double M,
                          void * params) {
  struct halos_workspace * params_local;
  struct background * pba;
  struct thermo * pth;
  double mu_1, mu_2, z, rho_c_over_rho_m;
  double one_halo_boost,factor_p_wave, F_of_c;
  params_local = params;
  pba = params_local->pba;
  pth = params_local->pth;
  z = params_local->z;
  mu_1 = log(1.+c)-c*pow(1.+c,-1);
  rho_c_over_rho_m = 1. + pba->Omega0_lambda/((pba->Omega0_cdm+pba->Omega0_b)*pow(1.+z,3.));
  if (pth->UCMH_DM_ann_type == s_wave) {
    mu_2 = (1./3.)*(1.-pow(1.+c,-3));
    one_halo_boost = (pth->Delta_c*rho_c_over_rho_m/3.)*pow(c,3)*mu_2*pow(mu_1,-2);
  } else {
    factor_p_wave = pow(_G_*H_per_H0(z,pba)*pba->h*1.e2*_km_over_Mpc_*M*_Sun_mass_over_kg_/pow(_c_,3),2./3.);
    factor_p_wave *= 4.*_PI_*pow(4.*pow(pth->Delta_c,4.),1./3.)*rho_c_over_rho_m;
    factor_p_wave *= pow(c,4.)/pow(mu_1,3.);
    F_of_c = -li2(-c)-6.*li3(1./(1.+c))-6.*log(1.+c)*li2(1./(1.+c))+6.*_zeta3_-(53./6.)
    +pow(log(1.+c),2)*(3.*log(c/(1.+c))+((2.+3.*c)/(c*(1.+c)))-0.5*(1.+pow(c,-2.)))
    +log(1.+c)*(pow(1.+c,-2.)+((1.+7.*c)/(c*(1.+c))))+pow(1.+c,-1)*(7.+pow(1.+c,-1)+(1./3.)*pow(1.+c,-2));
    one_halo_boost = factor_p_wave*F_of_c; // CHECK
  }
 return one_halo_boost;
}

/* define concentration function using prescription à la Macciò et al. (arXiv:0805.1926)*/
double c_NFW(double M,
             void * params) {
  /* M needs to be passed in units of solar masses */
struct halos_workspace * params_local;
struct background * pba;
double k_200=3.9, zF, z;
params_local = params;
pba = params_local->pba;
zF  = zF_wo_spike(params_local, 0.01*M);
z = params_local->z;
return k_200*pow(H_per_H0(zF,pba)/H_per_H0(z,pba),2./3.);
}

/* define Hubble parameter (divided by its present value H0) in terms of z, for the LCDM model  */
double H_per_H0(double z,
                struct background * pba) {
return sqrt(pba->Omega0_lambda+(pba->Omega0_cdm+pba->Omega0_b)*pow(1.+z,3)+(pba->Omega0_g+pba->Omega0_ur)*pow(1.+z,4));
}

/* computes the effective redshift for the formation of a halo of mass M  */
double zF_wo_spike(void * params, double M) {
 struct halos_workspace * params_local;
 struct background * pba;
 struct precision  * ppr;
 double rho_m_0,S, zF_no_spike;
 double gamma=6.*pow(_PI_,2);
 params_local = params;
 pba = params_local->pba;
 ppr = params_local->ppr;
 rho_m_0 = (pba->Omega0_cdm+pba->Omega0_b)*2.775e11*pow(pba->h,2); /* present matter density in units of M_sol/Mpc^3 */
 params_local->R = pow(M/(gamma*rho_m_0),1./3.); /* R is in units of Mpc, so k should be in units of Mpc^{-1}  */
 params_local->S0 = integrate_simpson(ppr->k_min, 1./params_local->R,ppr->Number_k, LOG, integrand_for_S,params_local);
 zF_no_spike = -1.+sqrt(params_local->S0)/_delta_crit_; //Here it is assumed for simplicity that D_growth(z)=1./(1.+z)
 return zF_no_spike;
}

double halo_function_high_mass(void * params,
                               double M) { /* M should be passed in units of solar mass, M_sol */
 struct halos_workspace * params_local;
 struct precision  * ppr;
 struct background * pba;
 struct thermo * pth;
 double z, R,omega_of_z, S_tot,rho_m_0;
 double nu, fnu, dS0_dM, P_R_standard,cut_off=1.;
 double H100=1./3000.; /* #To express Hubble constant in units of 100 Mpc^{-1} (having set c=1) */
 double gamma=6.*pow(_PI_,2);
 params_local = params;
 ppr = params_local->ppr;
 pba = params_local->pba;
 pth = params_local->pth;
 z = params_local->z;
 rho_m_0 = (pba->Omega0_cdm+pba->Omega0_b)*2.775e11*pow(pba->h,2); /* present matter density in units of M_sol/Mpc^3 */
 params_local->R = pow(M/(gamma*rho_m_0),1./3.); /* R is in units of Mpc, so k should be in units of Mpc^{-1}  */
 omega_of_z = _delta_crit_*D_growth(0.,pba)/D_growth(z,pba);
 S_tot = integrate_simpson(ppr->k_min, 1./params_local->R ,ppr->Number_k, LOG, integrand_for_S,params_local);
 nu = omega_of_z/sqrt(S_tot);
 fnu = sqrt(2./_PI_)*exp(-pow(nu,2)/2.);
 P_R_standard = pth->A_s*pow(cut_off,2)*pow(1./(pth->k_pivot*params_local->R),pth->n_s-1.);
 dS0_dM = (1./(3.*M))*(4./25.)*pow(pow(1./params_local->R,2.)*D_growth(0.,pba)*Transfer_f(1./params_local->R, pba,pth)/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h*pow(H100,2.)),2.)*P_R_standard;
 return (1./M)*(nu/(2.*S_tot))*fnu*dS0_dM;
}

/* define GROWTH FACTOR for the LCDM model */
double D_growth(double z,
                struct background * pba) {
 double Om_Lambda, Om_m;
 Om_m = (pba->Omega0_cdm+pba->Omega0_b)*pow(1.+z,3)/(pba->Omega0_lambda+(pba->Omega0_cdm+pba->Omega0_b)*pow(1.+z,3));
 Om_Lambda = pba->Omega0_lambda/(pba->Omega0_lambda+(pba->Omega0_cdm+pba->Omega0_b)*pow(1.+z,3));
 return pow(1.+z,-1)*(5.*Om_m/2.)/(pow(Om_m,4./7.) - Om_Lambda + (1.+Om_m/2.)*(1.+Om_Lambda/70.)); //Taken from formula A4 in arXiv: 9709112
}

/* compute concentration function for UCMHs */
double c_UCMH(double M,
              void * params) {
/* M needs to be passed in units of solar masses */
struct halos_workspace * params_local;
struct background * pba;
struct thermo * pth;
double conc_UCMH;
double z, zF, Omega_m_at_zF;
params_local = params;
pba = params_local->pba;
pth = params_local->pth;
z = params_local->z;
if (pth->M_at_Mthres == YES) {
  zF  = zF_Mthres_avg(params_local,z);
} else {
  zF  = zF_w_spike(params_local, M);
}
Omega_m_at_zF = (pba->Omega0_cdm+pba->Omega0_b)*pow(1.+zF,3.)/((pba->Omega0_cdm+pba->Omega0_b)*pow(1.+zF,3.) + pba->Omega0_lambda);
params_local->c3_over_mu = (3.*pth->f_2*Omega_m_at_zF/pth->Delta_c)*pow(H_per_H0(zF,pba)/H_per_H0(z,pba),2.);
conc_UCMH = find_root_Newton(pow(params_local->c3_over_mu,1./3.),1.e-2,g_UCMH,g_UCMH_prime, params_local);
return conc_UCMH;
}

/* compute concentration function for UCMHs */
double c_UCMH2(double zF,
              void * params) {
/* M needs to be passed in units of solar masses */
struct halos_workspace * params_local;
struct background * pba;
struct thermo * pth;
double conc_UCMH;
double z, Omega_m_at_zF;
params_local = params;
pba = params_local->pba;
pth = params_local->pth;
z = params_local->z;
Omega_m_at_zF = (pba->Omega0_cdm+pba->Omega0_b)*pow(1.+zF,3.)/((pba->Omega0_cdm+pba->Omega0_b)*pow(1.+zF,3.) + pba->Omega0_lambda);
params_local->c3_over_mu = (3.*pth->f_2*Omega_m_at_zF/pth->Delta_c)*pow(H_per_H0(zF,pba)/H_per_H0(z,pba),2.);
conc_UCMH = find_root_Newton(pow(params_local->c3_over_mu,1./3.),1.e-2,g_UCMH,g_UCMH_prime, params_local);
return conc_UCMH;
}

double g_UCMH(void * params,
              double c) {
  struct halos_workspace * params_local;
  params_local = params;
  double mu;
  mu = 2.*asinh(sqrt(c))-2.*sqrt(c/(1.+c));
 return pow(c,3) -params_local->c3_over_mu*mu;
}

double g_UCMH_prime(void * params,
                    double c) {
  struct halos_workspace * params_local;
  params_local = params;
  double mu_deriv;
  mu_deriv = c/sqrt(c*pow(c+1.,3));
 return 3.*pow(c,2) -params_local->c3_over_mu*mu_deriv;
}

double zF_w_spike(void * params,
                  double M) {
 struct halos_workspace * params_local;
 struct precision  * ppr;
 struct background * pba;
 struct thermo * pth;
 double rho_m_0,S_tot, S_standard, S_spike, zF;
 double H100=1./3000.; /* #To express Hubble constant in units of 100 Mpc^{-1} (having set c=1) */
 double gamma=6.*pow(_PI_,2);
 params_local = params;
 ppr = params_local->ppr;
 pba = params_local->pba;
 pth = params_local->pth;
 rho_m_0 = (pba->Omega0_cdm+pba->Omega0_b)*2.775e11*pow(pba->h,2); /* present matter density in units of M_sol/Mpc^3 */
 params_local->R = pow(M/(gamma*rho_m_0),1./3.); /* R is in units of Mpc, so k should be in units of Mpc^{-1}  */
 S_spike = (4./25.)*pth->A_spike*pow(pow(pth->k_spike,2)*D_growth(0.,pba)*Transfer_f(pth->k_spike, pba,pth)/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h*pow(H100,2)),2);
 S_standard = integrate_simpson(ppr->k_min, 1./params_local->R, ppr->Number_k, LOG, integrand_for_S,params_local);
 S_tot = S_spike + S_standard;
 if (pth->consider_only_spike_UCMH == _TRUE_) { //Here it is assumed for simplicity that D_growth(z)=1/(1+z)
   zF = -1.+sqrt(S_spike)/_delta_crit_;
 } else {
   zF = -1.+sqrt(S_tot)/_delta_crit_;
 }
 return zF;

}


double zF_Mthres_avg(void * params,
                     double z) {
 struct halos_workspace * params_local;
 struct precision  * ppr;
 struct background * pba;
 struct thermo * pth;
 double rho_m_0,S_tot, S_standard, S_spike;
 double nu_minus,nu_plus,omega_of_z, N_frac_spike, z_avg;
 double z_aprox;
 double H100=1./3000.; /* #To express Hubble constant in units of 100 Mpc^{-1} (having set c=1) */
 params_local = params;
 ppr = params_local->ppr;
 pba = params_local->pba;
 pth = params_local->pth;
 S_spike =  (4./25.)*pth->A_spike*pow(pow(pth->k_spike,2.)*D_growth(0.,pba)*Transfer_f(pth->k_spike, pba,pth)/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h*pow(H100,2.)),2.);
 S_standard = integrate_simpson(ppr->k_min, pth->k_spike ,ppr->Number_k, LOG, integrand_for_S, params_local);
 S_tot = S_spike + S_standard;
 omega_of_z = _delta_crit_*D_growth(0.,pba)/D_growth(z,pba);
 if (pth->consider_only_spike_UCMH == _TRUE_) {
   nu_minus = omega_of_z/sqrt(S_spike);
   nu_plus = 1.e10; //CHECK
   N_frac_spike = erfc(nu_minus/sqrt(2.));
 } else {
   nu_minus = omega_of_z/sqrt(S_tot);
   nu_plus  = omega_of_z/sqrt(S_standard);
   N_frac_spike = erfc(nu_minus/sqrt(2.)) - erfc(nu_plus/sqrt(2.));
 }
 z_avg = integrate_simpson(nu_minus, nu_plus,100000, LOG, integrand_for_zf_avg, params_local);
 z_avg = z_avg/N_frac_spike;
 if (N_frac_spike <1.e-100) { //to avoid divergences when N_frac_spike is tiny, substitute full expression of z_avg for approximate expression coming from a Taylor expansion (at nu_minus, nu_plus >> 1)
 z_avg =  -1. + ((1.+z)/nu_minus)*(1.-pow(nu_minus/nu_plus,2)*exp(-(pow(nu_plus,2)-pow(nu_minus,2))/2.))/(1.-(nu_minus/nu_plus)*exp(-(pow(nu_plus,2)-pow(nu_minus,2))/2.));
 }
 class_test(isnan(z_avg) || isinf(z_avg),
            pth->error_message,
            "The average formation redshift z_avg diverges,this should never happen");
return z_avg;
}

double integrand_for_zf_avg(void * params,
                            double nu) {
  struct halos_workspace * params_local;
  struct background * pba;
  double fnu, z, integrand;
  params_local = params;
  pba = params_local->pba;
  z = params_local->z;
  fnu = sqrt(2./_PI_)*exp(-pow(nu,2)/2.);
//  integrand =(pow(nu*D_growth(z,pba),-1)-1.)*fnu;
  integrand =(pow(nu/(1.+z),-1)-1.)*fnu; //Here it is assumed that D(z) = 1./(1.+z) for simplicity
  return integrand;
}

double integrand_boost_spike(void * params,
                             double zF) {
  struct halos_workspace * params_local;
  struct background * pba;
  struct thermo * pth;
  double nu,h_nu, c, integrand;
  double S_spike, S_standard, S_tot;
  double H100=1./3000.; /* #To express Hubble constant in units of 100 Mpc^{-1} (having set c=1) */
  params_local = params;
  pba = params_local->pba;
  pth = params_local->pth;
  S_spike =  (4./25.)*pth->A_spike*pow(pow(pth->k_spike,2.)*Transfer_f(pth->k_spike, pba,pth)/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h*pow(H100,2.)),2.);
  nu = _delta_crit_*(1.+zF)/sqrt(S_spike); //Here it is assumed that the growth factor is D(z)=1/(1+z)
  h_nu = pow(2.*_PI_,-2)*pow(3.,-3./2.)*nu*exp(-nu*nu/2.)*f_BBKS(nu);
  c = c_UCMH2(zF, params_local);
  integrand = pow(1.+zF,-1)*h_nu*one_halo_boost_UCMH2(c, zF, params_local);
  return integrand;
}


/* define function of concentration appearing for UCMH profiles */
double one_halo_boost_UCMH(double c,
                           double M,
                           void * params) {
  struct halos_workspace * params_local;
  struct background * pba;
  struct thermo * pth;
  double D, z, zF, Delta_t, rho_m_0;
  double mu_1, mu_2, rho_c_over_rho_m;
  double m_WIMP, sigmav_WIMP, one_halo_boost;
  double factor_p_wave, F_of_c;
  params_local = params;
  pba = params_local->pba;
  pth = params_local->pth;
  z = params_local->z;
  if (pth->M_at_Mthres == YES) {
    zF  = zF_Mthres_avg(params_local,z);
  } else {
    zF  = zF_w_spike(params_local, M);
  }
  if (z < zF) {
    Delta_t = MAX(tH0(z, pba)-tH0(zF, pba),tH0(zF, pba)/2.);
  } else {
    Delta_t = tH0(zF, pba)/2.;
  }
  Delta_t *= pow(100.*pba->h*1.02e-3,-1.); //Note: 1km/s/Mpc = 1.02e-3 Gyr^-1
  // following conversion is required, because following formulas assume m_WIMP in units of 1 TeV and sigmav in units of 3x10^{-26} cm^3/s
  m_WIMP = pth->DM_mass/1000.;
  sigmav_WIMP = pth->annihilation_cross_section/3.0e-26;
  D = pow(6.72e-11,-1)*pow(1.+zF,-2)*pow((m_WIMP/sigmav_WIMP)*(30./pth->f_2)*(13.8/Delta_t),2./3.);
  class_test((pow(D,-1.)> c),
             pth->error_message,
             "D^{-1} =%e is bigger than c =%e ,this should never happen",1./D, c);
  mu_1 = 2.*asinh(sqrt(c))-2.*sqrt(c/(1.+c));
  rho_c_over_rho_m = 1. + pba->Omega0_lambda/((pba->Omega0_cdm+pba->Omega0_b)*pow(1.+z,3.));

  if (pth->UCMH_DM_ann_type == s_wave) {
    mu_2 = (1./3.)+(2.*c+3.)/(2.*pow(1.+c,2))+log(c/(1.+c))-(2.*pow(D,-1)+3.)/(2.*pow(1.+pow(D,-1),2))+log(1.+D);
    one_halo_boost = (pth->Delta_c*rho_c_over_rho_m/3.)*pow(c,3)*mu_2*pow(mu_1,-2);
  } else {
    factor_p_wave = pow(_G_*H_per_H0(z,pba)*pba->h*1.e2*_km_over_Mpc_*M*_Sun_mass_over_kg_/pow(_c_,3),2./3.);
    factor_p_wave *= 4.*_PI_*pow(4.*pow(pth->Delta_c,4.),1./3.)*rho_c_over_rho_m;
    factor_p_wave *= pow(c,4.)/pow(mu_1,3.);
    F_of_c = -(8./15.)*pow(c,-5./2.)*pow(1.+c,-3./2.)
    *(-c*(-3.+c*(9.+c*(87.+70.*c)))+7.*pow(_PI_,2.)*(sqrt(pow(c,5.)*(1.+c))+sqrt(pow(c,7.)*(1.+c)))
    +3.*asinh(sqrt(c))*(-2.*sqrt(c*(1.+c))+5.*sqrt(pow(c,3.)*(1.+c))+12.*sqrt(pow(c,5.)*(1.+c))
    +(1.-16.*sqrt(pow(c,5.)*(1.+c))-16.*sqrt(pow(c,7.)*(1.+c))+c*(1.+2.*c)*(-1.+8.*c*(1.+c)))*asinh(sqrt(c))
    +8.*pow(c,5./2.)*pow(1.+c,3./2.)*(log(1.-exp(-2.*asinh(sqrt(c))))-5.*log(1.+exp(-2.*asinh(sqrt(c))))))
    +6.*(sqrt(pow(c,5.)*(1.+c))+sqrt(pow(c,7.)*(1.+c)))*(5.*li2(exp(-4.*asinh(sqrt(c)))) -12.*li2(exp(-2.*asinh(sqrt(c))))));
    one_halo_boost = factor_p_wave*F_of_c; //CHECK
  }

 return one_halo_boost;
}


/* define function of concentration appearing for UCMH profiles */
double one_halo_boost_UCMH2(double c,
                            double zF,
                            void * params) {
  struct halos_workspace * params_local;
  struct background * pba;
  struct thermo * pth;
  double D, z, Delta_t;
  double mmu_2, m_WIMP, sigmav_WIMP;
  params_local = params;
  pba = params_local->pba;
  pth = params_local->pth;
  z = params_local->z;
//  Delta_t = MAX(tH0(z, pba)-tH0(zF, pba),tH0(zF, pba)/2.);
// I checked that the full formula above gives essentially same results as approximate expression below
  Delta_t = tH0(z, pba);
  Delta_t *= pow(100.*pba->h*1.02e-3,-1.); //Note: 1km/s/Mpc = 1.02e-3 Gyr^-1
  // following conversion is required, because following formulas assume m_WIMP in units of 1 TeV and sigmav in units of 3x10^{-26} cm^3/s
  m_WIMP = pth->DM_mass/1000.;
  sigmav_WIMP = pth->annihilation_cross_section/3.0e-26;
  D = pow(6.72e-11,-1)*pow(1.+zF,-2)*pow((m_WIMP/sigmav_WIMP)*(30./pth->f_2)*(13.8/Delta_t),2./3.);
  class_test((pow(D,-1.)> c),
             pth->error_message,
             "D^{-1} =%e is bigger than c =%e ,this should never happen",1./D, c);
//  mmu_2 = (2./3.)+(2.*c+3.)/(pow(1.+c,2))+2.*log(c/(1.+c))-(2.*pow(D,-1)+3.)/(pow(1.+pow(D,-1),2))+2.*log(1.+D);
// I checked that the full formula above (that takes into account full dependence with c) gives essentially same results as approximate expression below
  mmu_2 = (2./3.)-3.+2.*log(D);
 return 2.*_PI_*pow(17.,2)*pow(1.+zF,3)*mmu_2;
}


/* define cosmic time (multiplied by H0) in terms of z, for the LCDM model (radiation and curvature not taken into account) */
double tH0(double z,
           struct background * pba) {
return (2./3.)*pow(pba->Omega0_lambda, -1./2.)*asinh(pow(pba->Omega0_lambda/(pba->Omega0_cdm+pba->Omega0_b), 1./2.)*pow(1.+z,-3./2.));
}


double integrand_boost_low_mass(void * params,
                                double M) {
  struct halos_workspace * params_local;
  double c, integrand;
  params_local = params;
  c = c_UCMH(M, params_local);
  integrand = M*halo_function_low_mass(params_local, M)*one_halo_boost_UCMH(c, M, params_local);
  return integrand;
}

double halo_function_low_mass(void * params,
                              double M) {
  struct halos_workspace * params_local;
  struct precision  * ppr;
  struct background * pba;
  struct thermo * pth;
  double z, R,omega_of_z, S_tot, S_standard, S_spike, rho_m_0;
  double nu, fnu, dS0_dM, P_R_standard, cut_off=1.;
 double H100=1./3000.; /* #To express Hubble constant in units of 100 Mpc^{-1} (having set c=1) */
 double gamma=6.*pow(_PI_,2);
 params_local = params;
 ppr = params_local->ppr;
 pba = params_local->pba;
 pth = params_local->pth;
 z = params_local->z;
 rho_m_0 = (pba->Omega0_cdm+pba->Omega0_b)*2.775e11*pow(pba->h,2); /* present matter density in units of M_sol/Mpc^3 */
 params_local->R = pow(M/(gamma*rho_m_0),1./3.);
 omega_of_z = _delta_crit_*D_growth(0.,pba)/D_growth(z,pba);
 S_spike = (4./25.)*pth->A_spike*pow(pow(pth->k_spike,2.)*D_growth(0.,pba)*Transfer_f(pth->k_spike, pba,pth)/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h*pow(H100,2.)),2.);
 S_standard = integrate_simpson(ppr->k_min, 1./params_local->R ,ppr->Number_k, LOG, integrand_for_S,params_local);
 S_tot = S_spike + S_standard;
 nu = omega_of_z/sqrt(S_tot);
 fnu = sqrt(2./_PI_)*exp(-pow(nu,2)/2.);
 P_R_standard = pth->A_s*pow(cut_off,2)*pow(1./(pth->k_pivot*params_local->R),pth->n_s-1.);
 dS0_dM = (1./(3.*M))*(4./25.)*pow(pow(1./params_local->R,2.)*D_growth(0.,pba)*Transfer_f(1./params_local->R, pba,pth)/((pba->Omega0_cdm+pba->Omega0_b)*pba->h*pba->h*pow(H100,2.)),2.)*P_R_standard;
 return (1./M)*(nu/(2.*S_tot))*fnu*dS0_dM;
}


/* function defined in eq. A.15 of J. M. Bardeen, J. R. Bond, N. Kaiser, and A. S. Szalay, Astrophys. J. 304, 15 (1986).*/
double f_BBKS(double x) {
return (pow(x,3)-3.*x)*(erf(sqrt(5./2.)*x)+erf(sqrt(5./2.)*x/2.))/2.
        +sqrt(2./(5.*_PI_))*(((31.*x*x/4.)+(8./5.))*exp(-5.*x*x/8.)+ ((x*x/2.)-(8./5.))*exp(-5.*x*x/2.));
}



double Transfer_f(double k,
                  struct background * pba,
                  struct thermo * pth) {
double kfs,Transfer,suppression;

if (pth->add_baryons_UCMH == _TRUE_) {
  Transfer = T_Hu(k, pba);
} else {
  Transfer = T_Hu_no_baryon(k, pba);
}

if (pth->add_suppression_kfs_UCMH == _TRUE_) {
  suppression = (1.-(2./3.)*pow(k/pth->k_fs,2))*exp(-pow(k/pth->k_fs,2)); //eq 47 in arXiv:0503387
} else {
  suppression = 1.;
}

Transfer *= suppression;
return Transfer;
}


double li2(double x) {

   const double PI = 3.1415926535897932;
   const double P[] = {
      0.9999999999999999502e+0,
     -2.6883926818565423430e+0,
      2.6477222699473109692e+0,
     -1.1538559607887416355e+0,
      2.0886077795020607837e-1,
     -1.0859777134152463084e-2
   };
   const double Q[] = {
      1.0000000000000000000e+0,
     -2.9383926818565635485e+0,
      3.2712093293018635389e+0,
     -1.7076702173954289421e+0,
      4.1596017228400603836e-1,
     -3.9801343754084482956e-2,
      8.2743668974466659035e-4
   };

   double y = 0, r = 0, s = 1;

   /* transform to [0, 1/2] */
   if (x < -1) {
      const double l = log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5*l - log(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const double l = log1p(-x);
      y = x/(x - 1);
      r = -0.5*l*l;
      s = -1;
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - log(x)*log(y);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const double l = log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(log(y) + 0.5*l);
      s = 1;
   } else {
      const double l = log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
   }

   const double y2 = y*y;
   const double y4 = y2*y2;
   const double p = P[0] + y * P[1] + y2 * (P[2] + y * P[3]) +
                    y4 * (P[4] + y * P[5]);
   const double q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3]) +
                    y4 * (Q[4] + y * Q[5] + y2 * Q[6]);

   return r + s*y*p/q;
}

/// Li_3(x) for x in [-1,0]
static double li3_neg(double x) {
   const double cp[] = {
      0.9999999999999999795e+0, -2.0281801754117129576e+0,
      1.4364029887561718540e+0, -4.2240680435713030268e-1,
      4.7296746450884096877e-2, -1.3453536579918419568e-3
   };
   const double cq[] = {
      1.0000000000000000000e+0, -2.1531801754117049035e+0,
      1.6685134736461140517e+0, -5.6684857464584544310e-1,
      8.1999463370623961084e-2, -4.0756048502924149389e-3,
      3.4316398489103212699e-5
   };

   const double x2 = x*x;
   const double x4 = x2*x2;
   const double p = cp[0] + x*cp[1] + x2*(cp[2] + x*cp[3]) +
      x4*(cp[4] + x*cp[5]);
   const double q = cq[0] + x*cq[1] + x2*(cq[2] + x*cq[3]) +
      x4*(cq[4] + x*cq[5] + x2*cq[6]);

   return x*p/q;
}


/// Li_3(x) for x in [0,1/2]
static double li3_pos(double x) {
   const double cp[] = {
      0.9999999999999999893e+0, -2.5224717303769789628e+0,
      2.3204919140887894133e+0, -9.3980973288965037869e-1,
      1.5728950200990509052e-1, -7.5485193983677071129e-3
   };
   const double cq[] = {
      1.0000000000000000000e+0, -2.6474717303769836244e+0,
      2.6143888433492184741e+0, -1.1841788297857667038e+0,
      2.4184938524793651120e-1, -1.8220900115898156346e-2,
      2.4927971540017376759e-4
   };

   const double x2 = x*x;
   const double x4 = x2*x2;
   const double p = cp[0] + x*cp[1] + x2*(cp[2] + x*cp[3]) +
      x4*(cp[4] + x*cp[5]);
   const double q = cq[0] + x*cq[1] + x2*(cq[2] + x*cq[3]) +
      x4*(cq[4] + x*cq[5] + x2*cq[6]);

   return x*p/q;
}

/**
 * @brief Real trilogarithm \f$\operatorname{Li}_3(x)\f$
 * @param x real argument
 * @return \f$\operatorname{Li}_3(x)\f$
 * @author Alexander Voigt
 */
double li3(double x) {
   const double zeta2 = 1.6449340668482264;
   const double zeta3 = 1.2020569031595943;

   // transformation to [-1,0] and [0,1/2]
   if (x < -1) {
      const double l = log(-x);
      return li3_neg(1/x) - l*(zeta2 + 1.0/6*l*l);
   } else if (x == -1) {
      return -0.75*zeta3;
   } else if (x < 0) {
      return li3_neg(x);
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5) {
      return li3_pos(x);
   } else if (x == 0.5) {
      return 0.53721319360804020;
   } else if (x < 1) {
      const double l = log(x);
      return -li3_neg(1 - 1/x) - li3_pos(1 - x)
         + zeta3 + l*(zeta2 + l*(-0.5*log(1 - x) + 1.0/6*l));
   } else if (x == 1) {
      return zeta3;
   } else if (x < 2) {
      const double l = log(x);
      return -li3_neg(1 - x) - li3_pos(1 - 1/x)
         + zeta3 + l*(zeta2 + l*(-0.5*log(x - 1) + 1.0/6*l));
   } else { // x >= 2.0
      const double l = log(x);
      return li3_pos(1/x) + l*(2*zeta2 - 1.0/6*l*l);
   }
}
