//Sampler for the mediation model with multiple (continuous) indicators for M and/or Y 
data {
  int n; //number of observations
  int M_ind;
  int Y_ind;
  vector[2] A_M; //prior standard deviation of first regression coefficients
  vector[3] A_Y; //prior standard deviation of second regression coefficients
  vector[n] X;
  vector[n] m_star[M_ind];  // a matrix with dimensions (M_ind * n) , that is M_ind vectors of length n
  vector[n] y_star[Y_ind];
}
parameters {
  real beta_1; real beta_2;
  /////////////////////
  real beta_3;
  /////////////////////
  real beta_0_M; real beta_0_Y;
  vector[2] lambda[M_ind-1];
  vector[2]  tau[Y_ind-1];
  vector<lower=0>[M_ind] ssq_m_star;
  vector<lower=0>[Y_ind] ssq_y_star;
  real<lower=0> ssq_M; real<lower=0> ssq_Y; 
  vector[n] M; vector [n] Y;

}
model {
  
  // Prior
  beta_0_M ~ normal(0,A_M[1]);
  beta_0_Y ~ normal(0,A_Y[1]);
  beta_1 ~ normal(0,A_M[2]);
  beta_2 ~ normal(0,A_Y[2]);
  /////////////////////
  beta_3 ~ normal(0,A_Y[3]);
  /////////////////////
  // M ~ normal(0,.000001); Y ~ normal(0,.000001);   //////CHECK THIS!!!
  if (M_ind > 1){
  for (i in 1:(M_ind-1))
  lambda[i] ~ normal(0,10);
  }
  if (Y_ind > 1){
  for (i in 1:(Y_ind-1))
  tau[i] ~ normal(0,10);
  }
  ssq_M ~ scaled_inv_chi_square(1, 3);
  ssq_Y ~ scaled_inv_chi_square(1, 3);
  ssq_m_star ~ scaled_inv_chi_square(1,3); //10000000, sqrt(2)); //scaled_inv_chi_square(10, 1);
  ssq_y_star ~ scaled_inv_chi_square(1,3); 
  
  // Likelihood
 M ~ normal(beta_0_M + beta_1*X, sqrt(ssq_M));
 /////////////////////
 Y ~ normal(beta_0_Y + beta_2*M + beta_3*X, sqrt(ssq_Y));
 // Y ~ normal(beta_0_Y + beta_2*M , sqrt(ssq_Y));
 /////////////////////  
 m_star[1] ~ normal(M, sqrt(ssq_m_star[1]));
 y_star[1] ~ normal(Y, sqrt(ssq_y_star[1]));
 if (M_ind > 1){
 for (i in 1:(M_ind-1))
  m_star[i+1] ~ normal(lambda[i,1] + lambda[i,2]*M, sqrt(ssq_m_star[i+1]));
 }
 if (Y_ind > 1){
 for (i in 1:(Y_ind-1))
  y_star[i+1] ~ normal(tau[i,1] + tau[i,2]*Y, sqrt(ssq_y_star[i+1])); 
 }
}
