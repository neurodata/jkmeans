using namespace arma;

class reducedDimARC {
 public:
  int n, d, p, K;
  mat Y;

  mat X;
  mat EV;
  mat EVV;
  mat Sigma;
  mat Zmu;
  double sigma2;
  bool randomStart;
  bool useEstep;

  double loglik;

  std::shared_ptr<Mixture> mix;
  unsigned int ver;

  reducedDimARC(mat _Y, int _d, int _K, bool _randomStart, double sigma2_ini,
                unsigned int _ver, bool fixW, bool _useEstep) {
    ver = _ver;

    Y = _Y;
    d = _d;
    K = _K;
    n = Y.n_rows;
    p = Y.n_cols;
    randomStart = _randomStart;
    useEstep = _useEstep;

    {
      mat U;
      vec s;
      mat V;
      svd(U, s, V, Y);
      U = U.cols(span(0, d - 1));
      s = s(span(0, d - 1));
      X = U * diagmat(s);
    }

    if (randomStart) {
      X += randn(n, d);
      // X = randn(n, d);
    }

    Sigma = eye(d, d);
    sigma2 = 1;

    mix.reset(new Mixture(X, K));

    // use kmeans++ for meansInput

    mat meansInput = kmeanspp(X, K);

    mix->initialize(meansInput, false, fixW, sigma2_ini);

    mix->Expectation();
    mix->Maximization();

    Zmu = (mix->zeta) * (mix->mu);
  };

  void Estep() {
    mat XX = X.t() * X;
    mat XY = X.t() * Y;
    mat XXinv = inv(XX + eye(d, d) * 1E-5);

    // expectation or maximize
    EV = XXinv * XY;

    // expectation
    if (useEstep) EVV = p * sigma2 * XXinv + EV * EV.t();
    // maximize
    else
      EVV = EV * EV.t();
  }

  void Mstep() {
    // sigma2
    sigma2 = (accu(Y % Y) - 2 * accu(vectorise(X.t()) % vectorise(EV * Y.t())) +
              accu(vectorise(X.t()) % vectorise(EVV * X.t()))) /
             n / p;
    // compute X
    mat invSigma = inv(Sigma);
    mat m = Y * EV.t() / sigma2 + Zmu * invSigma;
    mat v = inv(EVV / sigma2 + invSigma);
    X = m * v;

    switch (ver) {
      case 1:
        mix->runEM(100, 1E-5);
      case 2: {
        mix->Expectation();
        mix->Maximization();
      }
    }

    Zmu = (mix->zeta) * (mix->mu);
    Sigma = diagmat(mix->Sigma);
  }

  double computeLoglik() {
    mat XZmu = X - Zmu;
    double logdetSigma;
    double sign;
    log_det(logdetSigma, sign, Sigma);

    return accu(XZmu % trans(solve(Sigma, XZmu.t()))) + n * logdetSigma +
           n * p * log(sigma2) + n * p;
  }

  void runEM(int steps) {
    double cur_loglik = INFINITY;
    for (int i = 0; i < steps; ++i) {
      Estep();
      Mstep();
      double new_loglik = computeLoglik();
      cout << new_loglik << endl;
      if (new_loglik > cur_loglik & i > 1)
        break;
      else
        cur_loglik = new_loglik;
    }
    loglik = cur_loglik;
  }
};