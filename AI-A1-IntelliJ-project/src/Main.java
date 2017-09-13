public class Main {
  public static void main(String[] args) {
    InputHandler ih = new InputHandler();
    Model model = ih.readInitialModel();
    int[] obsSeq = ih.readObsSeq();
    int n = model.a.length;
    int m = model.b[0].length;
    BaumWelch bw = new BaumWelch(model, obsSeq, n, m);
    bw.run();
    MatrixUtils.prettyPrint(bw.model.a);
  }
}
