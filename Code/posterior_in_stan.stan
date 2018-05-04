data { 
  int N;
  int p;
  matrix[N, p] x;
  int y[N];
}

parameters {
  vector[p] beta;
}

model{
  // prior
  beta ~ normal(0, 100);
  
  //likelihood
  y ~ bernoulli_logit(x * beta);
}