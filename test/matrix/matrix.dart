library matrixtest;

import 'package:dolt/dolt.dart';
import 'package:unittest/unittest.dart';

part 'matrix1d.dart';
part 'matrix2d.dart';

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

  group('1d Selected Dense Matrix', () {
    var t = new SelectedDenseMatrix1dTests();
    t.testAll();
  });

  group('2d Dense Matrix', () {
    var t = new DenseMatrix2dTests();
    t.testAll();
  });

  group('2d Selected Dense Matrix', () {
    var t = new DenseMatrix2dTests();
    t.testAll();
  });

}
