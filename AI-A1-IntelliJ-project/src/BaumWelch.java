class Model {

}

public class BaumWelch {
  private Model model;
  private int[] obsSeq;
  private int n, m, t;
  private int[][] alpha;

  public BaumWelch(Model initialModel, int[] obsSeq, int n, int m) {
    this.model = initialModel;
    this.n = n;
    this.m = m;
    this.t = obsSeq.length;
    this.obsSeq = obsSeq;
  }

  private int[][] alpha_pass() {
    return null;
  }
}
