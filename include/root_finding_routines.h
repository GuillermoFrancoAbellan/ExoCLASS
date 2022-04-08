/* some basic root-finding routines */

double find_root_Newton(double x_guess,
                        double tolerance,
                        double (*func)(void * params_for_function, double x),
                        double (*func_prime)(void * params_for_function,double x),
                        void * params_for_function){
int i=0, N_max=100;
double x_new, x_old;
double f, f_prime;

x_old = x_guess;
f = (*func)(params_for_function, x_old);
f_prime = (*func_prime)(params_for_function, x_old);

while (fabs(f)>tolerance) {
  x_new = x_old - f/f_prime;
  x_old = x_new;
  f = (*func)(params_for_function, x_old);
  f_prime = (*func_prime)(params_for_function, x_old);
  i++;
  if (i>N_max) {
    printf("ERROR: Newton method to find root is not converging. \n");
    exit(EXIT_FAILURE);
  }

  if (f_prime == 0.0 ) {
    printf("ERROR: Newton method is crashing because f'(x_i)=0 \n");
    exit(EXIT_FAILURE);
  }

}

return x_old;
}
