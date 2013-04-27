part of matrixtest;

abstract class Matrix1dTests {
  Matrix1d fromList(List<double> values);

  void testAll() {
    test('dot product tests', testDot);
    test('aggregate method tests', testAggregate);
    test('assign method tests', testAssign);
  }

  void testDot() {
    var a = fromList([1.0, 2.0, 3.0, 5.0]);
    var b = fromList([1.0, 0.0, 3.0, 0.0]);

    expect(a.dot(b), equals(10.0));
  }

  void testAggregate() {
    var a = fromList([1.0, 2.0, 3.0, 5.0]);
    var sum = a.aggregate((cumulative, current) => cumulative + current, (x) => 2 * x);
    expect(sum, equals(22.0));

    var b = fromList([2.0, 3.0, 0.0, 4.0]);
    sum = a.aggregateOther(b, (cumulative, current) => cumulative + current, (x, y) => x * y);
    expect(sum, equals(28.0));
  }

  void testAssign() {
    var a = fromList([1.0, 1.0, 1.0]);
    a.assign(2.0);
    for (int i = 0; i < 3; i++) {
      expect(a[i], equals(2.0));
    }

    a.assign([1.0, 2.0, 3.0]);
    expect(a[0], equals(1.0));
    expect(a[1], equals(2.0));
    expect(a[2], equals(3.0));
  }
}

class SparseMatrix1dTests extends Matrix1dTests {
  Matrix1d fromList(List<double> values) {
    return new SparseMatrix1d(values);
  }
}

class DenseMatrix1dTests extends Matrix1dTests {
  Matrix1d fromList(List<double> values) {
    return new DenseMatrix1d(values);
  }
}

class SelectedSparseMatrix1dTests extends Matrix1dTests {
  Matrix1d fromList(List<double> values) {
    // Add some extra cruft to the values we get, and then create a view
    // around them.  It's a little cheesy, but suffices for a first pass.
    var allValues = [0.0, 42.0, 7.0]
        ..addAll(values)
        ..addAll([11.0, -67.3]);
    var offsets = [];
    var emap = new Map<int, double>();

    for (int i = 0; i < allValues.length; i++) {
      emap[i] = allValues[i];
      offsets.add(i);
    }
    return new SelectedSparseMatrix1d(values.length, emap, 0, 1, offsets, 3);
  }
}

class SelectedDenseMatrix1dTests extends Matrix1dTests {
  Matrix1d fromList(List<double> values) {
    // Add some extra cruft to the values we get, and then create a view
    // around them.  It's a little cheesy, but suffices for a first pass.
    var allValues = [0.0, 42.0, 7.0]
        ..addAll(values)
        ..addAll([11.0, -67.3]);
    var offsets = [];
    for (int i = 0; i < allValues.length; i++) {
      offsets.add(i);
    }
    return new SelectedDenseMatrix1d(values.length, allValues, 0, 1, offsets, 3);
  }
}
