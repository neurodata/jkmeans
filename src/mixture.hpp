using namespace arma;

class Mixture {
 public:
  mat y;
  int n;
  int p;
  int K;  // max components
  int J;  // truncation

  mat mu;
  double sigma2;
  vec w;
  mat zeta;

  vec gradient;
  mat hessian;

  Mixture(mat _y, int _K, int _J) {
    y = _y;
    K = _K;
    J = _J;
    n = y.n_rows;
    p = y.n_cols;

    // initialzie weights
    w = ones(K) / (double)K;
    sigma2 = 1;
    mu = ones<mat>(K, p);
    mu.randn();

    mat means;
    bool status = kmeans(means, trans(y), K, random_subset, 10, false);

    mu = means.t();

    zeta = ones<mat>(n, K);

    gradient = zeros(p * K);
    hessian = zeros<mat>(p * K, p * K);
  }

  double loglik(rowvec x, rowvec mu, double sigma2) {
    rowvec diff = x - mu;

    return -dot(diff, diff) / 2 / sigma2 - p / 2 * log(sigma2);
  }

  double compTotalLoglik() {
    vec marginal_loglik = zeros<vec>(n);

    for (int k = 0; k < K; ++k) {
      mat diff = (y - repmat(mu.row(k), n, 1));

      vec diff2 =
          exp(-sum(diff % diff, 1) / sigma2 / 2.0 - p / 2.0 * log(sigma2));

      marginal_loglik += w(k) * diff2;
    }

    return accu(log(marginal_loglik));
  }

  void Expectation() {
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < K; ++k) {
        rowvec y_local = y.row(i);
        rowvec mu_local = mu.row(k);
        zeta(i, k) = log(w(k)) + loglik(y_local, mu_local, sigma2);
        if (std::isnan(zeta(i, k))) {
          cout << w(k) << endl;
          cout << loglik(y_local, mu_local, sigma2) << endl;
          cout << y_local << endl;
          cout << mu_local << endl;
          cout << sigma2 << endl;
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

      local_zeta(find(local_zeta < local_zeta(trunc_idx))).fill(0);

      zeta.row(i) = local_zeta / accu(local_zeta);
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
      sum_total_weights += total_weights * p;

      // sigma2 = (accu(sum(diff % diff, 1) % weights)) /
      // (total_weights * p + 0.1);

      // compute w
      // w(k) = accu(weights) / (double)n;
    }

    sigma2 = (diff2 / sum_total_weights);
  }

  vec computeGradient(mat mu) {
    vec gradient(p * K);

    for (int k = 0; k < K; ++k) {
      mat diff = y - repmat(mu.row(k), n, 1);

      vec zeta_local = zeta.col(k);
      vec gradient_local = diff.t() * zeta_local / sigma2;

      double hessian_diag_local =
          -accu(zeta_local) / sigma2;  // same for all p elements

      int idx = k * p;
      int end_idx = idx + p - 1;
      gradient(span(idx, end_idx)) = gradient_local;
    }

    return gradient;
  }

  mat computeHessian() {
    mat hessian = zeros<mat>(p * K, p * K);

    for (int k = 0; k < K; ++k) {
      mat diff = y - repmat(mu.row(k), n, 1);

      vec zeta_local = zeta.col(k);

      double hessian_diag_local =
          -accu(zeta_local) / sigma2;  // same for all p elements

      int idx = k * p;
      int end_idx = idx + p - 1;

      for (int i = 0; i < p; ++i) {
        hessian(idx + i, idx + i) = hessian_diag_local;
      }
    }
    return hessian;
  }

  void updateSigma2() {
    double diff2 = 0;
    double sum_total_weights = 0;
    for (int k = 0; k < K; ++k) {
      mat diff = y - repmat(mu.row(k), n, 1);

      vec zeta_local = zeta.col(k);

      // lazy: use MLE to estimate sigma2 :P
      diff2 += accu(sum(diff % diff, 1) % zeta_local);
      double total_weights = accu(zeta_local);
      sum_total_weights += total_weights * p;
    }

    sigma2 = (diff2 / sum_total_weights);
  }

  void GradientDescent(int steps) {
    mat B = zeros<mat>(p * K, p * K);

    for (int i = 0; i < steps; ++i) {
      /* code */

      Expectation();

      double curLoglik = compTotalLoglik();

      mat hessian = computeHessian();

      mat old_mu = mu;


      vec Q_n_minus_1_bar_n_minus_1 = computeGradient(old_mu);

      // run one descent
      mat H = hessian - B;
      double factor = 0.5;
      while (rcond(H) < 1E-8) {
        H = hessian - factor * B;
        factor *= 0.5;
      }

      mu = old_mu - trans(reshape(solve(H, Q_n_minus_1_bar_n_minus_1), p, K));

      Expectation();

      vec Q_n_minus_1_bar_n = computeGradient(old_mu);

      vec g_n = Q_n_minus_1_bar_n - Q_n_minus_1_bar_n_minus_1;
      vec s_n = reshape(mu - old_mu, p * K, 1);

      vec v_n = g_n - B * s_n;
      double c_n = 1 / dot(v_n, s_n);
      B = B + c_n * v_n * v_n.t();

      updateSigma2();

      double newLoglik = compTotalLoglik();

      factor = 0.5;
      mat direction = mu - old_mu;
      while (newLoglik + 1E-8 < curLoglik) {
        mu = old_mu + factor * direction;
        factor *= 0.5;
        updateSigma2();
        newLoglik = compTotalLoglik();
        if (factor < 0.1) {
          mu = old_mu + 0.1 * direction;
          updateSigma2();
          break;
        }
      }

      if(std::isnan(newLoglik) || std::isinf(newLoglik) || mu.has_nan()){
        mu = old_mu;
        updateSigma2();
      }




      //cout << newLoglik << endl;
    }
  }

  void runEM(int steps, double tol) {

    mat mu0 = mu;

    for (int i = 0; i < steps; ++i) {
      Expectation();
      Maximization();

      mat diff = mu-mu0;

      if( accu(diff % diff) < tol)
        break;
      else
        mu0 = mu;

     // cout << compTotalLoglik() << endl;
    }
  }

  void runQNEM(int steps) {
    Expectation();

    GradientDescent(steps);

  }
};
