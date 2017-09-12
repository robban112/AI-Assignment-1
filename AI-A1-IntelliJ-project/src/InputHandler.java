/**
 * Created by RobertLorentz on 2017-09-12.
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class InputHandler {
  private Scanner scanner;

  InputHandler() {
    this.scanner = new Scanner(System.in);
  }

  public Model readInitialModel() {
    double[][] a = createMatrix(scanner.nextLine());
    double[][] b = createMatrix(scanner.nextLine());
    double[] pi = createMatrix(scanner.nextLine())[0];
//        int[] obs = createObsSeq(scanner.nextLine());
    return new Model(a, b, pi);
  }

  public int[] readObsSeq() {
    String line = scanner.nextLine();
    String[] obsStrs = line.split(" ");
    int obsNum = Integer.parseInt(obsStrs[0]);
    int[] obs = new int[obsNum];
    for(int i = 0; i<obsNum; i++) {
      obs[i] = Integer.parseInt(obsStrs[i+1]);
    }
    return obs;
  }

  private static double[][] createMatrix(String line) {
    String[] nums = line.split(" ");
    ArrayList<String> tail = new ArrayList<>();
    tail.addAll(Arrays.asList(nums));
    tail.remove(0);
    tail.remove(0);
    double[][] matrix = createMatrixFromList(nums[0], nums[1], tail);
    return matrix;
  }

  private static double[][] createMatrixFromList(String strRows, String strCols, ArrayList<String> elems) {
    int rows = Integer.parseInt(strRows);
    int cols = Integer.parseInt(strCols);
    double[][] arr = new double[rows][cols];
    int rc = 0;
    int cc = 0;
    for (String elem : elems) {
      if (cc >= cols) {
        cc=0;
        rc++;
      }
      arr[rc][cc] = Double.parseDouble(elem);
      cc++;
    }
    return arr;
  }
}
