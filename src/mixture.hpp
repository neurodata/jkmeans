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
  mat zeta;

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

    mat means;
    bool status = kmeans(means, trans(y), K, random_subset, 10, false);

    mu = means.t();

    zeta = ones<mat>(n, K);
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
        zeta(i, k) =
            log(w(k)) + loglik(y_local, mu_local, sigma2(k));
        if (std::isnan(zeta(i, k))) {
          cout << w(k) << endl;
          cout << loglik(y_local, mu_local, sigma2(k)) << endl;
          cout << y_local << endl;
          cout << mu_local << endl;
          cout << sigma2(k) << endl;
        }
      }
      rowvec local_zeta = zeta.row(i);

      local_zeta -= local_zeta.max();

      // if(local_zeta.has_nan()){
      //   cout<< local_zeta<<endl;
      //   cout<< zeta<<endl;
      // }

      uvec indices = sort_index(local_zeta, "descend");

      int trunc_idx = indices(J - 1);

      local_zeta = exp(local_zeta);

      local_zeta(
          find(local_zeta < local_zeta(trunc_idx)))
          .fill(0);

      zeta.row(i) =
          local_zeta / accu(local_zeta);
    }
  }

  void Maximization() {

    double diff2 = 0;
    double sum_total_weights = 0;

    for (int k = 0; k < K; ++k) {
      vec weights = zeta.col(k);
      double total_weights = accu(weights);

      // compute mu
      if (total_weights > 0)
        mu.row(k) = sum(y % repmat(weights, 1, p), 0) / total_weights;
      else
        mu.row(k) = trans(randn(p));

      // compute sigma2
      mat diff = (y - repmat(mu.row(k), n, 1));

      diff2 += accu(sum(diff % diff, 1) % weights);
      sum_total_weights += total_weights * p ;

      // sigma2(k) = (accu(sum(diff % diff, 1) % weights)) /
                  // (total_weights * p + 0.1);

      // compute w
      // w(k) = accu(weights) / (double)n;
    }

    sigma2.fill(diff2 / sum_total_weights);
  }

  void run(int steps) {
    for (int i = 0; i < steps; ++i) {
      Expectation();
      Maximization();
    }
  }
};