using namespace arma;

class Mixture {
 public:
  mat y;
  int n;
  int p;
  int K;  // max components
  int J;  // truncation

  mat mu;
  vec sigma2;
  vec w;
  mat posterior_weight;

  Mixture(mat _y, int _K, int _J) {
    y = _y;
    K = _K;
    J = _J;
    n = y.n_rows;
    p = y.n_cols;

    // initialzie weights
    w = ones(K) / (double)K;
    sigma2 = ones(K);
    mu = ones<mat>(K, p);
    mu.randn();

    posterior_weight = ones<mat>(n, K);
  }

  double loglik(rowvec x, rowvec mu, double sigma2) {
    rowvec diff = x - mu;

    return -dot(diff, diff) / 2 / sigma2 - p / 2 * log(sigma2);
  }

  void Expectation() {
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < K; ++k) {
        rowvec y_local = y.row(i);
        rowvec mu_local = mu.row(k);
        posterior_weight(i, k) =
            log(w(k)) + loglik(y_local, mu_local, sigma2(k));
        if (std::isnan(posterior_weight(i, k))) {
          cout << w(k) << endl;
          cout << loglik(y_local, mu_local, sigma2(k)) << endl;
          cout << y_local << endl;
          cout << mu_local << endl;
          cout << sigma2(k) << endl;
        }
      }
      rowvec local_posterior_weight = posterior_weight.row(i);

      local_posterior_weight -= local_posterior_weight.max();

      // if(local_posterior_weight.has_nan()){
      //   cout<< local_posterior_weight<<endl;
      //   cout<< posterior_weight<<endl;
      // }

      uvec indices = sort_index(local_posterior_weight, "descend");

      int trunc_idx = indices(J - 1);

      local_posterior_weight = exp(local_posterior_weight);

      local_posterior_weight(
          find(local_posterior_weight < local_posterior_weight(trunc_idx)))
          .fill(0);

      posterior_weight.row(i) =
          local_posterior_weight / accu(local_posterior_weight);
    }
  }

  void Maximization() {
    for (int k = 0; k < K; ++k) {
      vec weights = posterior_weight.col(k);
      double total_weights = accu(weights);

      // compute mu
      if (total_weights > 0)
        mu.row(k) = sum(y % repmat(weights, 1, p), 0) / total_weights;
      else
        mu.row(k) = trans(randn(p));

      // compute sigma2
      mat diff = (y - repmat(mu.row(k), n, 1));

      sigma2(k) = (accu(sum(diff % diff, 1) % weights) + 0.1) /
                  (total_weights * p + 0.1);

      // compute w
      w(k) = accu(weights) / (double)n;
    }
  }

  void run(int steps) {
    for (int i = 0; i < steps; ++i) {
      Expectation();
      Maximization();
    }
  }
};