using namespace arma;

class Mixture {
 public:
  mat y;
  int n;
  int p;
  int K;  // max components
  int J;  // truncation

  mat mu;
  vec Sigma;
  vec w;
  mat zeta;

  bool fixW;
  bool normalizeZeta;
  bool flexJ;
  double zetaTrunc;

  Mixture(mat _y, int _K, int _J) {
    y = _y;
    K = _K;
    J = _J;
    n = y.n_rows;
    p = y.n_cols;

    // initialize weights
    w = ones(K) / (double)K;
    Sigma = ones(p);
    mu = ones<mat>(K, p);
    mu.randn();

    zeta = ones<mat>(n, K);
  }

  void initialize(mat meansInput, bool useKmeansIni, bool _fixW, bool _flexJ,
                  double _zetaTrunc, double sigma2_ini, bool _normalizeZeta) {
    Sigma.fill(sigma2_ini);

    if (useKmeansIni) {
      mat means;
      bool status = kmeans(means, trans(y), K, static_spread, 10, false);
      mu = means.t();

      // initialize sigma2 via K-means
      int _J = J;  // back up the user defined J
      J = 1;
      Expectation();

      updateSigma();

      J = _J;

    } else {
      mu = meansInput;
    }

    normalizeZeta = _normalizeZeta;
    fixW = _fixW;
    flexJ = _flexJ;
    zetaTrunc = _zetaTrunc;
  }

  double loglik(rowvec x, rowvec mu, vec Sigma) {
    rowvec diff = x - mu;

    return -accu(diff % diff / Sigma.t()) / 2.0 - 1.0 / 2.0 * accu(log(Sigma));
  }

  double compTotalLoglik() {
    double marginal_loglik = 0;

    for (int k = 0; k < K; ++k) {
      mat diff = (y - repmat(mu.row(k), n, 1));

      mat diff2 = diff % diff;

      double l = exp(-accu(diff2.each_row() / Sigma.t()) / 2.0 -
                     1 / 2.0 * accu(log(Sigma)));

      marginal_loglik += w(k) * l;
    }

    return accu(log(marginal_loglik));
  }

  void Expectation() {
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < K; ++k) {
        rowvec y_local = y.row(i);
        rowvec mu_local = mu.row(k);
        zeta(i, k) = log(w(k)) + loglik(y_local, mu_local, Sigma);
        if (std::isnan(zeta(i, k))) {
          cout << w(k) << endl;
          cout << loglik(y_local, mu_local, Sigma) << endl;
          cout << y_local << endl;
          cout << mu_local << endl;
          // cout << Sigma << endl;
        }
      }
      rowvec local_zeta = zeta.row(i);

      local_zeta -= local_zeta.max();

      if (!flexJ) {
        uvec indices = sort_index(local_zeta, "descend");

        int trunc_idx = indices(J - 1);

        local_zeta = exp(local_zeta);

        local_zeta = local_zeta / accu(local_zeta);

        local_zeta(find(local_zeta < local_zeta(trunc_idx))).fill(0);
      }
      if (flexJ) {
        local_zeta = exp(local_zeta);
        local_zeta /= accu(local_zeta);
        local_zeta(find(local_zeta < zetaTrunc)).fill(0);
      }
      if (normalizeZeta)
        zeta.row(i) = local_zeta / accu(local_zeta);
      else
        zeta.row(i) = local_zeta;
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

      // // compute sigma2
      // mat diff = (y - repmat(mu.row(k), n, 1));

      // diff2 += accu(sum(diff % diff, 1) % weights);
      // sum_total_weights += total_weights * p;

      // compute w
      if (!fixW) w(k) = accu(weights) / (double)n;
    }

    updateSigma();
  }

  void updateSigma() {
    vec diff2 = zeros(p);
    double sum_total_weights = 0;
    for (int k = 0; k < K; ++k) {
      mat diff = y.each_row() - mu.row(k);

      vec zeta_local = zeta.col(k);

      diff2 += (zeta_local.t() * (diff % diff)).t();

      // diff2 += accu(sum(diff % diff, 1) % zeta_local);
      double total_weights = accu(zeta_local);
      sum_total_weights += total_weights;
    }

    Sigma = (diff2 / sum_total_weights);
  }

  void runEM(int steps, double tol) {
    mat mu0 = mu;

    double cur_loglik = -INFINITY;

    for (int i = 0; i < steps; ++i) {
      Expectation();
      Maximization();

      mat diff = mu - mu0;

      if (accu(diff % diff) < tol)
        break;
      else
        mu0 = mu;

      // cout << compTotalLoglik() << endl;
    }
  }

  void sortBy1stDinMu() {
    uvec indices = sort_index(mu.col(0));
    zeta = zeta.cols(indices);
    mu = mu.rows(indices);
    w = w(indices);
  }

  uvec clusteringMAP() { return arma::index_max(zeta, 1); }
};
