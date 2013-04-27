library matrixtest;

import 'package:dolt/dolt.dart';
import 'package:unittest/unittest.dart';

part 'matrix1d.dart';

void matrixTestAll() {
  group('1d Sparse Matrix', () {
    var t = new SparseMatrix1dTests();
    t.testAll();
  });

  group('1d Selected Sparse Matrix', () {
    var t = new SelectedSparseMatrix1dTests();
    t.testAll();
  });

  group('1d Dense Matrix', () {
    var t = new DenseMatrix1dTests();
    t.testAll();
  });
}
