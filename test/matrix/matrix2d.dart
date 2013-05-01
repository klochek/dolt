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