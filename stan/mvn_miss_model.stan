data {
    int N; // number of participants
    int J; // number of observations per participant
    int K; // number of treatment groups
    matrix[N, J] y; // outcomes
    array[N] int x; // treatment indicator
    array[N, J] int R;  // missingness indicator
    matrix[K, K-1] Xdes; // design
    real med;
    real mad;
}

transformed data {
  int<lower=0> n_missing = 0;
  for(n in 1:N) {
      for(j in 1:J) {
          n_missing += R[n, j] == 1;
      }
  }
}

parameters {
    real alpha;
    vector[J-1] b_visit;
    matrix[K-1,J-1] b_treat;
    cholesky_factor_corr[J] L_Omega;
    vector<lower=0>[J] sigma;
    array[n_missing] real y_missing;
}

transformed parameters {
    cholesky_factor_cov[J] L_Sigma;
    corr_matrix[J] Omega;
    cov_matrix[J] Sigma;
    // cell means by visit and treatment group
    matrix[K, J] mu;
    matrix[K, J-1] beta;
    for(j in 1:J) {
        for(k in 1:K) {
            if(j == 1) {
                mu[k, j] = alpha;
            } else {
                beta[k,j-1] = Xdes[k] * b_treat[:, j-1];
                mu[k,j] = alpha + b_visit[j-1] + beta[k,j-1];
            }
        }
    }
    // shared covariance of visits across all treatment groups
    L_Sigma = diag_pre_multiply(sigma, L_Omega);
    Omega = multiply_lower_tri_self_transpose(L_Omega);
    Sigma = quad_form_diag(Omega, sigma);
}

model { 
    int pos_miss = 1;
    matrix[N, J] y_imp;

    // priors
    L_Omega ~ lkj_corr_cholesky(2);
    sigma ~ student_t(3, 0, 5);
    alpha ~ normal(med, mad);
    b_visit ~ normal(0, 5);
    to_vector(b_treat) ~ normal(0, 5);

    // likelihood
    for(n in 1:N) {
        for(j in 1:J) {
            if (R[n, j] == 1) {
                y_imp[n, j] = y_missing[pos_miss];
                pos_miss += 1;
            } else {
                y_imp[n, j] = y[n, j];
            }
        }
        y_imp[n] ~ multi_normal(mu[x[n]], Sigma);
    }
}
