/* some basic integration routines */

double integrate_simpson_nonuniform(double * x_array,
                                    double * y_array,
                                    int first_index, //by default will be 0
                                    int Num_samples) { //Number of elements of x_array and y_array to integrate
                                                     // it can be smaller than the length of these arrays, in case we want
                                                     // to integrate only some interval, that's why first_index is useful

double integral=0.;
double h0, h1, hph, hdh, hmh;
int Num_intervals, i;
double * h_step;

Num_intervals = Num_samples-1;
h_step = malloc(Num_intervals*sizeof(double));
for (i=first_index; i <first_index+Num_intervals; i++) {
  h_step[i] = x_array[i+1]-x_array[i];
}

for (i=first_index+1; i < first_index+Num_intervals; i =i+2) {
  h0 = h_step[i-1];
  h1 = h_step[i];
  hph = h1+h0;
  hdh = h1/h0;
  hmh = h1*h0;
  integral += (hph/6.)*((2.-hdh)*y_array[i-1]+(hph*hph/hmh)*y_array[i]+(2.-1./hdh)*y_array[i+1]);
}

if (Num_intervals%2 == 1) {
  h0 = h_step[first_index+Num_intervals-2];
  h1 = h_step[first_index+Num_intervals-1];
  integral += y_array[first_index+Num_intervals]*(2.*h1*h1+3.*h0*h1)/(6.*(h0+h1));
  integral += y_array[first_index+Num_intervals-1]*(h1*h1+3.*h1*h0)/(6.*h0);
  integral -= y_array[first_index+Num_intervals-2]*h1*h1*h1/(6.*h0*(h0+h1));
}

free(h_step);

return integral;
}
