/* some basic integration routines */

enum integration_sampling {LIN, LOG};


double integrate_simpson(double xmin,
                         double xmax,
                         int N,
                         enum integration_sampling sampling,
                         double (*function)(void * params_for_function, double x),
                         void * params_for_function) {
double integral, h;
double sum_odd = 0., sum_even = 0.;
double * x;
double * f;
int i;

x = malloc((N+1)*sizeof(double));
f = malloc((N+1)*sizeof(double));

if ((xmin <= 0.) && (sampling == LOG)) {
 printf("ERROR: you cannot use negative or zero values with the LOG sampling \n");
 exit(EXIT_FAILURE);
}

if (sampling == LIN) {
  h = (xmax-xmin)/N;
} else {
  h = (log10(xmax)-log10(xmin))/N;
}

for (i=0; i <= N; i++) {
 if (sampling == LIN) {
   x[i] = xmin+i*h;
   f[i] = (*function)(params_for_function, x[i]);
 } else {
   x[i] = xmin*pow(10,i*h);
   f[i] = (*function)(params_for_function, x[i])*x[i]*log(10); // since dx = dln(x)*x = dlog_10(x)*ln(10)*x
 }
}

for (i=1; i < N; i++) {
  if (i%2 == 1) {
    sum_odd += f[i];
  } else {
    sum_even += f[i];
  }
}

integral = (h/3.)*(f[0]+4.*sum_odd+ 2.*sum_even+f[N]);
free(x);
free(f);

return integral;
}
