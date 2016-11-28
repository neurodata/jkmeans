using namespace arma;

mat kmeanspp(mat X, int K) {
  int n = X.n_rows;
  uvec C(K);

  std::random_device rd;
  std::mt19937 gen(rd());

  vec pr = ones(n);

  auto randDiscrete = [&](vec pr) {
    std::vector<double> p_stdvec = conv_to<std::vector<double> >::from(pr);
    std::discrete_distribution<int> distribution(p_stdvec.begin(),
                                                 p_stdvec.end());

    return distribution(gen);
  };

  C(0) = randDiscrete(pr);

  mat D = X.each_row() - X.row(C(0));
  vec D2 = sum(D % D, 1);

  for (int i = 1; i < K; i++) {
    pr = D2;
    pr(C(span(0, i - 1))).fill(0);
    C(i) = randDiscrete(pr);

    mat D = X.each_row() - X.row(C(i));
    vec D2new = sum(D % D, 1);

    D2 = min(D2, D2new);
  }
  return X.rows(C);
}