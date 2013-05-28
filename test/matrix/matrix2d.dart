part of matrixtest;

abstract class Matrix2dTests {
  Matrix2d fromList(List<List<double>> values);

  void testAll() {
    test('2dx2d multiplication tests', testMultiplication);
  }

  void testMultiplication() {
    var a = fromList([[3.0, 2.0], [7.0, 1.0]]);
    var b = fromList([[2.0, 5.0], [4.0, 3.0]]);

    var c = a.mul2d(b, null);

    expect(c.get(0, 0), equals(14.0));
    expect(c.get(1, 1), equals(38.0));
    expect(c.get(1, 0), equals(18.0));
    expect(c.get(0, 1), equals(21.0));
  }
}

class DenseMatrix2dTests extends Matrix2dTests {
  Matrix2d fromList(List<List<double>> values) {
    var r = new DenseMatrix2d(values);

    // Set a tiny block size, so that we can stress the blocking part of multiplication.
    r.blockSize = 1;
    return r;
  }
}

class SelectedDenseMatrix2dTests extends Matrix2dTests {
  Matrix2d fromList(List<List<double>> values) {
    var flatVals = [];
    for (int i = 0; i < values.length; i++) {
      for (int j = 0; j < values[i].length; j++) {
        flatVals.add(values[i][j]);
      }
    }
    var r = new SelectedDenseMatrix2d(flatVals, [0, 2], [0, 1], 0);
    return r;
  }
}

class SparseMatrix2dTests extends Matrix2dTests {
  Matrix2d fromList(List<List<double>> values) {
    var r = new SparseMatrix2d(values);
    return r;
  }
}

class SelectedSparseMatrix2dTests extends Matrix2dTests {
  Matrix2d fromList(List<List<double>> values) {
    var flatVals = new Map<int, double>();
    for (int i = 0; i < values.length; i++) {
      for (int j = 0; j < values[i].length; j++) {
        flatVals[i * values[i].length + j] = values[i][j];
      }
    }
    var r = new SelectedSparseMatrix2d(flatVals, [0, 2], [0, 1], 0);
    return r;
  }
}
