public class Main {
  public static void main(String[] args) {
    InputHandler ih = new InputHandler();
    Model model = ih.readInitialModel();
    int[] obsSeq = ih.readObsSeq();
    int n = model.a.length;
    int m = model.b[0].length;
//    for (double[] row : model.a) {
//      System.out.println();
//      for (double val : row)
//        System.out.print(val + " ");
//    }
    BaumWelch bw = new BaumWelch(model, obsSeq, n, m);
    bw.run();
//    for (double val : MatrixUtils.getCol(bw.model.a, 1)) {
//      System.out.print(val + " ");
//    }
    for (double[] row : bw.model.a) {
      System.out.println();
      for (double val : row)
        System.out.print(val + " ");
    }
  }
}
