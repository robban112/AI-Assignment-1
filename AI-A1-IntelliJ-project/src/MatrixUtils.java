public class MatrixUtils {
  public static double[] getCol(double[][] m, int n) {
    double[] col = new double[m.length];
    for (int i = 0; i < m.length; i++) {
      col[i] = m[i][n];
    }
    return col;
  }

  public static double[] vectorMult(double[] v1, double[] v2) {
    assert v1.length == v2.length;

    double[] resVec = new double[v1.length];
    for (int i = 0; i < v1.length; i++) {
      resVec[i] = v1[i]*v2[i];
    }
    return resVec;
  }
}
