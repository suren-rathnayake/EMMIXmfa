void cgmf(double *X, double *A, double *B, int *samps, int* vars,
          int *q,  double *lambda, int* maxit,
          double* Z0, double *cor_rate) //, double *verr, double *grr)
{
int  k = q[0];
double gamma = lambda[0];
double error = 0;
double Err = 0;
double alpha, S, beta, Z;
int f, l, i, j, p, n, n_iti = maxit[0];

beta = cor_rate[0]; 
Z    = Z0[0]; 
n    = samps[0]; 
p    = vars[0];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (l = 1; l <= n_iti; l = l + 1) { // global iterations
  for (i = 0; i <= p - 1; i = i + 1) { // itirate over variables
    for (j = 0; j <= n - 1; j = j + 1) { // itirate over samples

      //x = X[i + p * j];
      S = 0.0;

      for (f = 0; f <= k - 1; f = f + 1) {
        S = S + A[i + p * f] * B[j + n * f];
      }

      error = X[i + p * j] - S;
      for (f = 0; f <= k - 1; f = f + 1) {
        alpha        =  A[i + p * f] * B[j + n * f];
        A[i + p * f] +=  gamma * error * B[j + n * f];
        error        +=  alpha - A[i + p * f] * B[j + n * f];
        alpha        =  A[i + p * f] * B[j + n * f];
        B[j + n * f] +=  gamma * error * A[i + p * f];
        error        +=  alpha - A[i + p * f] * B[j + n * f];
      }
    }
  }

  Err = 0.0;
  S = 0.0;
  for (i = 0; i <= p-1; i = i + 1) {
    for (j = 0; j <= n-1; j = j + 1) {
      for (f = 0; f <= k-1; f = f + 1) {
        S = S + A[i + p * f] * B[j + n * f];
      }

      Err = Err + (X[i + p * j] - S) * (X[i + p * j] - S);
      S = 0.0;
    }
  }

  Err = Err / (n * p);
  //verr[l - 1] = Err;
  //grr[l - 1] = gamma;
  if (Err > Z) {

    gamma = beta * gamma;
  } else {

    Z = Err;
  }
}

//Z0[0] = Err;
}

