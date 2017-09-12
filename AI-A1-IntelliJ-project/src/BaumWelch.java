public class BaumWelch {
  private Model model;
  private int[] obsSeq;
  private int n, m, maxT;
  private double[][] alpha, beta;
  private double[][] gamma;
  private double[][][] digamma;
  private double[] cs;

  public BaumWelch(Model initialModel, int[] obsSeq, int n, int m) {
    this.model = initialModel;
    this.n = n;
    this.m = m;
    this.maxT = obsSeq.length;
    this.obsSeq = obsSeq;
    this.alpha = new double[maxT][n];
    this.beta = new double[maxT][n];
    this.digamma = new double[maxT-1][n][n];
    this.gamma = new double[maxT][n];
    this.cs = new double[maxT];
  }

  private void alphaPass() {
    // compute alpha[0]
    double[] firstObsProbs = MatrixUtils.getCol(model.b, obsSeq[0]);
    alpha[0] = MatrixUtils.vectorMult(model.pi, firstObsProbs);

    // empty prev scaling values
    for (int i = 0; i < cs.length; i++)
      cs[i] = 0;

    // scale the alpha[0]
    for (double a : alpha[0])
      cs[0] += a;
    cs[0] = 1/cs[0];
    for (int i = 0; i < n; i++)
      alpha[0][i] *= cs[0];

    // compute alpha[0..T-1]
    for (int t = 1; t < maxT; t++) {
      for (int i = 0; i < n; i++) {
        alpha[t][i] = 0;
        for (int j = 0; j < n; j++)
          alpha[t][i] += alpha[t-1][j]*model.a[j][i]*model.b[i][obsSeq[t]];
        cs[t] += alpha[t][i];
      }

      // scale alpha[t]
      cs[t] = 1/cs[t];
      for (int i = 0; i < n; i++)
        alpha[t][i] *= cs[t];
    }
  }

  private void betaPass() {
    // initialize beta[T-1]
    for (int i = 0; i < n; i++)
      beta[maxT-1][i] = cs[maxT-1];

    // compute beta[0..T-2]
    for (int t = maxT-2; t <= 0; t--) {
      for (int i = 0; i < n; i++) {
        beta[t][i] = 0;
        for (int j = 0; j < n; j++)
          beta[t][i] = beta[t+1][j]*model.a[i][j]*model.b[j][obsSeq[t+1]];
        // scale beta[t][i]
        beta[t][i] *= cs[t];
      }
    }
  }

  private void gammaDigamma() {
      for (int t = 0; t < maxT-1; t++)
        for (int i = 0; i < n; i++) {
          gamma[t][i] = 0;
          double denom = 0;
          for (int j = 0; j < n; j++)
            denom += alpha[t][i] * model.a[i][j] * beta[t+1][j];
          for (int j = 0; j < n; j++) {
            digamma[t][i][j] = (alpha[t][i]*model.a[i][j]*beta[t+1][j])/denom;
            gamma[t][i] += digamma[t][i][j];
          }
        }
  }

  private void reestimatePi() {
    model.pi = gamma[0];
  }

  private void reestimateA() {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double numer = 0;
        double denom = 0;
        for (int t = 0; t < maxT - 1; t++) {
          numer += digamma[t][i][j];
          denom += gamma[t][i];
        }
        model.a[i][j] = numer / denom;
      }
    }
  }

  private void reestimateModel() {
    gammaDigamma();
    reestimatePi();
    reestimateA();
  }

  private double computeLogProb() {
    double logProb = 0;
    for (double c : cs)
      logProb += Math.log(c);
    return logProb;
  }

  public void run() {
    int maxIters = 30;
    int iters = 0;
    double oldLogProb = -1000;
    double logProb = -999;
    while (iters < maxIters && logProb > oldLogProb) {
      oldLogProb = logProb;
      alphaPass();
      betaPass();
      reestimateModel();
      logProb = computeLogProb();
      iters++;
    }
  }

  private void reestimateB() {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        double numer = 0;
        double denom = 0;
        for (int t = 0; t < maxT; t++) {
          denom += gamma[t][i];
          if (obsSeq[t] == j) {
            numer += gamma[t][i];
          }
        }
        model.b[i][j] = numer / denom;
      }
    }
  }
}
