public class Main {
  public static void main(String[] args) {
    InputHandler ih = new InputHandler();
    Model model = ih.readInitialModel();
    int[] obsSeq = ih.readObsSeq();
    int n = model.a.length;
    int m = model.b[0].length;
    BaumWelch bw = new BaumWelch(model, obsSeq, n, m);
    bw.run();
    printModel(bw.model);
    //MatrixUtils.prettyPrint(bw.model.a);
  }

  public static void printModel(Model model) {
    printMatrix(model.a);
    printMatrix(model.b);
  }

  public static void printMatrix(double[][] matrix) {
    System.out.print(matrix.length + " " + matrix[0].length + " ");
    for (int r = 0; r<matrix.length; r++) {
      for (int c = 0; c<matrix[r].length; c++) {
        System.out.print(matrix[r][c]);
        if(r != matrix.length-1 || c != matrix[r].length-1) {
          System.out.print(" ");
        }
      }
    }
    System.out.println();
  }
}
