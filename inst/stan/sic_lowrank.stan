data {

  int<lower=1> N;
  int<lower=1> Q;
  int<lower=1> P;
  int<lower=1> R;
  int<lower=0, upper = 1> y_star; // input as c(y_star)
  matrix[N, P] y_lagged;
  matrix[N, Q] X;

  real beta_prior_mean;
  real<lower=0> beta_prior_sd;
  real UV_prior_mean;
  real<lower=0> UV_prior_sd;

}

parameters {

  matrix[P, R] U;
  matrix[P, R] V_I;
  matrix[P, R] V_C;
  matrix[Q, P] Beta_I;
  matrix[Q, P] Beta_C;

}

transformed parameters {

  matrix[N, P] Lambda;
  matrix[N, P] y_probs;

  Lambda =
  exp(
    (1.0 - y_lagged) .* ( X * Beta_I + y_lagged * (U * V_I') )  +
    (y_lagged) .* ( X * Beta_C + y_lagged * (U * V_C') )
  );
  y_probs =
  1.0 - exp(-Lambda);

}

model {

  // for (p in 1:P) {
  //   for (n in 1:N) {
  //     y_star[n,p] ~ bernoulli(y_probs[n,p]);
  //   }
  // }
  y_star ~ bernoulli(to_vector(y_probs));

  to_vector(Beta_I) ~ normal(beta_prior_mean,beta_prior_sd);
  to_vector(Beta_C) ~ normal(beta_prior_mean,beta_prior_sd);
  to_vector(U) ~ normal(UV_prior_mean,UV_prior_sd);
  to_vector(V_I) ~ normal(UV_prior_mean,UV_prior_sd);
  to_vector(V_C) ~ normal(UV_prior_mean,UV_prior_sd);

}

