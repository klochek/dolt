/*
 * Copyright � 1999 CERN - European Organization for Nuclear Research.
 * Permission to use, copy, modify, distribute and sell this software and its documentation for any
 * purpose is hereby granted without fee, provided that the above copyright notice appear in all
 * copies and that both that copyright notice and this permission notice appear in supporting
 * documentation. CERN makes no representations about the suitability of this software for any
 * purpose. It is provided "as is" without expressed or implied warranty.
 *
 * Dart port Copyright (c) 2013 Chris Klochek.
 */

part of dolt;

/**
 * Abstract base class for 2-d matrices.
 * A matrix has a number of cells, which are assigned upon instance construction.
 * Elements are accessed via zero based indexes
 */
abstract class Matrix2d {
  bool _isView = false;

  /**
   * Ensures that this matrix can hold at least the specified number of non-zero (non-null) cells
   * without needing to allocate new internal memory.  If necessary, allocates new internal memory
   * and increases the capacity of the matrix.
   */
  void _ensureCapacity(int minNonZeros) {
  }

  /**
   * Returns whether the receiver is a view or not.
   */
  bool get isView {
    return _isView;
  }

  /**
   * Releases any superfluous internal memory. An application can use this operation to minimize the
   * storage of the receiver.
   */
  void _trimToSize() {
  }


  /// The number of columns this matrix has.
  int _columns;

  /// The number of rows this matrix has.
  int _rows;

  /// The number of elements between two rows.
  int _rowStride;

  /// The number of elements between two columns.
  int _columnStride;

  /// The index of the first row element.
  int _rowZero;

  /// The index of the first column element.
  int _columnZero;


  /**
   * Returns the [int] position of the given absolute rank within the (virtual or non-virtual)
   * internal 1-dimensional array.
   */
  int _columnOffset(int absRank) {
    return absRank;
  }

  /**
   * Returns the [int] absolute rank of the given relative [rank].
   */
  int _columnRank(int rank) {
    return _columnZero + rank * _columnStride;
  }

  /**
   * Returns the [int] position of the given absolute rank within the (virtual or non-virtual)
   * internal 1-dimensional array.
   */
  int _rowOffset(int absRank) {
    return absRank;
  }

  /**
   * Returns the [int] absolute rank of the given relative [rank].
   */
  int _rowRank(int rank) {
    return _rowZero + rank * _rowStride;
  }

  /**
   * Checks whether the matrix contains the given box and throws a [RangeError], if necessary.
   */
  void _checkValidBox(int row, int column, int height, int width) {
    if (column < 0 || width < 0 || column + width > _columns ||
        row < 0 || height < 0 || row + height > _rows) {
      throw new RangeError('${toStringShort()}, column:$column, ' +
          'row:$row, width:$width, height:$height');
    }
  }

  /**
   * Sanity check for operations requiring a column index to be within bounds.
   */
  void _checkValidColumn(int column) {
    if (column < 0 || column >= _columns) {
      throw new RangeError('Attempted to access ${toStringShort()} at column=$column');
    }
  }

  /**
   * Checks whether indexes are legal and throws a [RangeError], if necessary.
   */
  void _checkValidColumnIndexes(List<int> indexes) {
    indexes.forEach(_checkValidColumn);
  }

  /**
   * Sanity check for operations requiring a row index to be within bounds.
   */
  void _checkValidRow(int row) {
    if (row < 0 || row >= _rows) {
      throw new RangeError('Attempted to access ${toStringShort()} at row=$row');
    }
  }

  /**
   * Checks whether indexes are legal and throws a [RangeError], if necessary.
   */
   void _checkValidRowIndexes(List<int> indexes) {
     indexes.forEach(_checkValidRow);
  }

  /**
   * Sanity check for operations requiring two matrices with the same number of columns and rows.
   */
  void _checkValidShape(Matrix2d b) {
    if (_columns != b._columns || _rows != b._rows) {
      throw new ArgumentError('Incompatible dimensions: ${toStringShort()}' +
          ' and ${b.toStringShort()}');
    }
  }

  /**
   * Sanity check for operations requiring matrices with the same number of columns and rows.
   */
  void _checkValidShapes(Matrix2d b, Matrix2d c) {
    if (_columns != b._columns || _rows != b._rows ||
        _columns != c._columns || _rows != c._rows) {
      throw new ArgumentError('Incompatible dimensions: ${toStringShort()}, ' +
          '${b.toStringShort()}, ${c.toStringShort()}');
    }
  }

  /**
   * Returns the number of columns.
   */
  int get columns {
    return _columns;
  }

  /**
   * Returns the number of rows.
   */
  int get rows {
    return _rows;
  }

  /**
   * Returns the position of the given coordinate within the (virtual or non-virtual) internal
   * 1-dimensional array.
   */
  int index(int row, int column) {
    return _rowOffset(_rowRank(row)) + _columnOffset(_columnRank(column));
  }

  /**
   * Sets up a matrix with a given number of [rows] and [columns].
   */
  void _setUp(int rows, int columns) {
    _setUpWithParams(rows, columns, 0, 0, columns, 1);
  }

  /**
   * Sets up a matrix with a given number of rows and columns and the given strides.
   */
  void _setUpWithParams(int rows, int columns, int rowZero, int columnZero, int rowStride,
      int columnStride) {
    if (rows < 0 || columns < 0) {
      throw new ArgumentError('negative size');
    }
    _rows = rows;
    _columns = columns;

    _rowZero = rowZero;
    _columnZero = columnZero;

    _rowStride = rowStride;
    _columnStride = columnStride;

    _isView = false;

    // TODO(klochek): Do we even want this, anymore?
    if (columns * rows > 4000000000) {
      throw new ArgumentError('matrix too large');
    }
  }

  /**
   * Returns the number of cells which is [rows * columns].
   */
  int get size {
    return _rows * _columns;
  }

  /**
   * Returns a string representation of the receiver's shape.
   */
  String toStringShort() {
    throw 'Not implemented yet!';
    //return AbstractFormatter.shape(this);
  }

  /**
   * Self modifying version of viewColumnFlip().
   */
  void vColumnFlip() {
    if (_columns > 0) {
      _columnZero += (_columns - 1) * _columnStride;
      _columnStride = -_columnStride;
      _isView = true;
    }
  }

  /**
   * Self modifying version of viewTranspose().
   */
  void vTranspose() {
    int tmp;

    tmp = _rows;
    _rows = _columns;
    _columns = tmp;

    tmp = _rowStride;
    _rowStride = _columnStride;
    _columnStride = tmp;

    tmp = _rowZero;
    _rowZero = _columnZero;
    _columnZero = tmp;

    _isView = true;
  }

  /**
   * Self modifying version of viewPart().
   */
  void vPart(int row, int column, int height, int width) {
    _checkValidBox(row, column, height, width);
    _rowZero += _rowStride * row;
    _columnZero += _columnStride * column;

    _rows = height;
    _columns = width;

    _isView = true;
  }

  /**
   * Self modifying version of viewRowFlip().
   */
  void vRowFlip() {
    if (_rows > 0) {
      _rowZero += (_rows - 1) * _rowStride;
      _rowStride = -_rowStride;
      _isView = true;
    }
  }

  /**
   * Self modifying version of viewStrides().
   */
  void vStrides(int rowStride, int columnStride) {
    if (rowStride <= 0 || columnStride <= 0) {
      throw new RangeError('illegal strides: $rowStride, $columnStride');
    }
    _rowStride *= rowStride;
    _columnStride *= columnStride;
    if (_rows != 0) {
      _rows = ((_rows - 1) / rowStride + 1).floor();
    }
    if (_columns != 0) {
      _columns = ((_columns - 1) / columnStride + 1).floor();
    }
    _isView = true;
  }


  double aggregate(aggr, f) {
    if (size == 0) {
      return double.NAN;
    }
    double a = f(get(rows - 1, columns - 1));
    int d = 1;
    for (int row = _rows - 1; row >= 0; row--) {
      for (int column = _columns - d - 1; column >= 0; column--) {
        a = aggr(a, f(get(row, column)));
      }
      d = 0;
    }
    return a;
  }

  double aggregateOther(Matrix2d other, aggr, f) {
    _checkValidShape(other);
    if (size == 0) {
      return double.NAN;
    }
    double a = f(get(rows - 1, columns - 1), other.get(rows - 1, columns - 1));
    int d = 1;
    for (int row = rows - 1; row >= 0; row--) {
      for (int column = columns - d - 1; column >= 0; column--) {
        a = aggr(a, f(get(row, column), other.get(row, column)));
      }
      d = 0;
    }
    return a;
  }

  /**
   * Sets all cells to the state specified by [values].
   */
   void assign(List<List<double>> values) {
    if (values.length != _rows) {
      throw new ArgumentError('Must have same number of rows: values=${values.length} rows=$_rows');
    }
    for (int row = rows - 1; row >= 0; row--) {
      List<double> currentRow = values[row];
      if (currentRow.length != columns) {
        throw new ArgumentError('Must have same number of columns in every row: ' +
            'currentRow=${currentRow.length} columns=$_columns');
      }
      for (int column = columns - 1; column >= 0; column--) {
        set(row, column, currentRow[column]);
      }
    }
  }

  /**
   * Sets all cells to the state specified by [value].
   */
  void assignValue(double value) {
    int r = _rows;
    int c = _columns;
    for (int row = 0; row < r; row++) {
      for (int column = 0; column < c; column++) {
        set(row, column, value);
      }
    }
  }

  /**
   * Assigns the result of a function to each cell.
   */
  void assignFunc(function) {
    for (int row = rows - 1; row >= 0; row--) {
      for (int column = columns - 1; column >= 0; column--) {
        set(row, column, function(get(row, column)));
      }
    }
  }

  /**
   * Replaces all cell values of the receiver with the values of another matrix.
   * Both matrices must have the same number of rows and columns.
   */
  void assignMatrix(Matrix2d other) {
    if (other == this) {
      return;
    }
    _checkValidShape(other);
    if (haveSharedCells(other)) {
      other = other.copy();
    }

    for (int row = rows - 1; row >= 0; row--) {
      for (int column = columns - 1; column >= 0; column--) {
        set(row, column, other.get(row, column));
      }
    }
  }

  /**
   * Assigns the result of a function to each cell.
   */
  void assignOtherFunc(Matrix2d y, function) {
    _checkValidShape(y);
    for (int row = rows - 1; row >= 0; row--) {
      for (int column = columns - 1; column >= 0; column--) {
        set(row, column, function(get(row, column), y.get(row, column)));
      }
    }
  }

  /**
   * Returns the number of cells having non-zero values; ignores tolerance.
   */
  int cardinality() {
    int cardinality = 0;
    for (int row = rows - 1; row >= 0; row--) {
      for (int column = columns - 1; column >= 0; column--) {
        if (get(row, column) != 0) {
          cardinality++;
        }
      }
    }
    return cardinality;
  }

  /**
   * Constructs and returns a deep copy of the receiver.
   */
   Matrix2d copy() {
    return like()
        ..assign(this);
  }

  /**
   * Returns whether all cells are equal to the given value.
   */
  bool equalsValue(double value) {
    throw 'Not implemented yet!';
    //return cern.colt.matrix.linalg.Property.DEFAULT.equals(this,value);
  }

  /**
   * Compares this object against the specified object.  Ignores tolerances.
   */
   bool equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (obj is! Matrix2d) {
      return false;
    }

    throw 'Not implemented yet!';
    //return cern.colt.matrix.linalg.Property.DEFAULT.equals(this,(DoubleMatrix2D) obj);
  }

  /**
   * Assigns the result of a function to each non-zero cell.
   */
  void forEachNonZero(function) {
    for (int row = rows - 1; row >= 0; row--) {
      for (int column = columns - 1; column >= 0; column--) {
        double value = get(row, column);
        if (value != 0) {
          double r = function(row, column, value);
          if (r != value) {
            set(row, column, r);
          }
        }
      }
    }
  }

  /**
   * Returns the matrix cell value at coordinate <tt>[row,column]</tt>.
   */
   double safeGet(int row, int column) {
     _checkValidBox(row, column, rows, columns);
    return get(row, column);
  }

  /**
   * Returns the content of this matrix if it is a wrapper; or [this] otherwise.
   */
   Matrix2d get content {
    return this;
  }

  /**
   * Fills the coordinates and values of cells having non-zero values into the specified lists.
   * Fills into the lists, starting at index 0.
   */
  void getNonZeros(List<int> rowList, List<int> columnList, List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int r = rows;
    int c = columns;
    for (int row = 0; row < r; row++) {
      for (int column = 0; column < c; column++) {
        double value = get(row, column);
        if (value != 0) {
          rowList.add(row);
          columnList.add(column);
          valueList.add(value);
        }
      }
    }
  }

  /**
   * Returns the [double] value in the matrix cell value at [row] and [column].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  double get(int row, int column);

  /**
   * Returns [true] if both matrices share at least one identical cell.
   */
  bool haveSharedCells(Matrix2d other) {
    if (other == null) {
      return false;
    }
    if (this == other) {
      return true;
    }
    return content.haveSharedCellsRaw(other.content);
  }

  /**
   * Returns [true] if both matrices share at least one identical cell.
   */
  bool haveSharedCellsRaw(Matrix2d other) {
    return false;
  }

  /**
   * Construct and returns a new empty matrix of the same dynamic type as this matrix, having the
   * same number of rows and columns.
   */
   Matrix2d like() {
    return likeWithDimensions(rows, columns);
  }

  /**
   * Construct and returns a new empty matrix of the same dynamic type as this matrix, having the
   * specified number of [rows] and [columns].
   */
  Matrix2d likeWithDimensions(int rows, int columns);

  /**
   * Construct and returns a new 1-d matrix of the corresponding dynamic type of this matrix, with
   * the specified [size].
   */
  Matrix1d like1d(int size);

  /**
   * Construct and returns a new 1-d matrix of the corresponding dynamic type of this matrix, with
   * specified [size], [zero], and [stride].
   */
  Matrix1d like1dWithDimensions(int size, int zero, int stride);

  /**
   * Sets the matrix cell at coordinate [row] and [column] to the specified [value].  Throws a
   * [RangeError] if no such index exists.
   */
  void safeSet(int row, int column, double value) {
    _checkValidBox(row, column, rows, columns);
    set(row, column, value);
  }

  /**
   * Sets the matrix cell at coordinate [row] [column] to the specified [value].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  void set(int row, int column, double value);

  /**
   * Constructs and returns a [List] of [List]s containing the cell values.  The outer-most list is
   * a list of rows.
   */
  List<List<double>> toArray() {
    var values = new List<List<double>>.fill(rows, new List<double>(columns));
    for (int row = rows - 1; row >= 0; row--) {
      var currentRow = values[row];
      for (int column = columns - 1; column >= 0; column--) {
        currentRow[column] = get(row, column);
      }
    }
    return values;
  }

  /**
   * Returns a string representation using default formatting.
   */
  String toString() {
    throw 'Not implemented!';
    //return new cern.colt.matrix.doublealgo.Formatter().toString(this);
  }

  /**
   * Constructs and returns a new view equal to the receiver.  The view is a shallow clone.
   */
  Matrix2d view() {
    throw 'Not implemented yet';
  }

  /**
   * Constructs and returns a new slice view representing the rows of the given column.  The
   * returned view is backed by this matrix, so changes in the returned view are reflected in this
   * matrix, and vice-versa.  To obtain a slice view on subranges, construct a sub-ranging view
   * via [viewPart], then apply this method to the sub-range view.
   */
  Matrix1d viewColumn(int column) {
    _checkValidColumn(column);
    int viewSize = _rows;
    int viewZero = index(0, column);
    int viewStride = _rowStride;
    return like1d(viewSize, viewZero, viewStride);
  }

  /**
   * Constructs and returns a new flip view along the column axis.  What used to be column 0 is now
   * column [columns - 1], and vice versa.
   * The returned view is backed by this matrix, so changes in the returned view are reflected in
   * this matrix, and vice-versa.
   */
  Matrix2d viewColumnFlip() {
    return view()
        ..vColumnFlip();
  }

  /**
   * Constructs and returns a new transposition view.  This is a zero-copy transposition, taking
   * constant time.  The returned view is backed by this matrix, so changes in the returned view are
   * reflected in this matrix, and vice-versa.
   */
  Matrix2d viewTranspose() {
    return view()
        ..vTranspose();
  }

  /**
   * Constructs and returns a new sub-range view that is a [height] x [width] sub matrix starting at
   * [row] and [column].  Operations on the returned view can only be applied to the restricted
   * range.  Any attempt to access coordinates not contained in the view will throw an [RangeError].
   */
  Matrix2d viewPart(int row, int column, int height, int width) {
    return view()
        ..vPart(row, column, height, width);
  }

  /**
   * Constructs and returns a new slice view representing the columns of the given row.  The
   * returned view is backed by this matrix, so changes in the returned view are reflected in this
   * matrix, and vice-versa.  To obtain a slice view on subranges, construct a sub-ranging view
   * via [viewPart], then apply this method to the sub-range view.
   */
  Matrix1d viewRow(int row) {
    _checkValidRow(row);
    int viewSize = _columns;
    int viewZero = index(row, 0);
    int viewStride = _columnStride;
    return like1dWithDimensions(viewSize, viewZero, viewStride);
  }

  /**
   * Constructs and returns a new flip view along the row axis.  What used to be row 0 is now row
   * [rows - 1], and vice versa.  The returned view is backed by this matrix, so changes in the
   * returned view are reflected in this matrix, and vice-versa.
   */
  Matrix2d viewRowFlip() {
    return view()
        ..vRowFlip();
  }

  /**
   * Constructs and returns a new selection view that is a matrix holding the indicated cells.
   * After construction, view.get(i, j) == this.get(rowIndexes[i], columnIndexes[j]).  Indexes can
   * occur multiple times and can be in arbitrary order.
   */
  Matrix2d viewSelection(List<int> rowIndexes, List<int> columnIndexes) {
    if (rowIndexes == null) {
      rowIndexes = new List<int>(rows);
      for (int i = rows - 1; i >= 0; i--) {
        rowIndexes[i] = i;
      }
    }
    if (columnIndexes == null) {
      columnIndexes = new List<int>(columns);
      for (int i = columns - 1; i >= 0; i--) {
        columnIndexes[i] = i;
      }
    }

    _checkValidRowIndexes(rowIndexes);
    _checkValidColumnIndexes(columnIndexes);
    var rowOffsets = new List<int>(rowIndexes.length);
    var columnOffsets = new List<int>(columnIndexes.length);
    for (int i = rowIndexes.length - 1; i >= 0; i--) {
      rowOffsets[i] = _rowOffset(_rowRank(rowIndexes[i]));
    }
    for (int i = columnIndexes.length - 1; i >= 0; i--) {
      columnOffsets[i] = _columnOffset(_columnRank(columnIndexes[i]));
    }
    return viewSelectionLike(rowOffsets, columnOffsets);
  }

  /**
   * Constructs and returns a new selection view that is a matrix holding all rows matching the
   * given condition.  To match columns, use a transpose view.
   */
  Matrix2d viewSelectionOnCondition(condition) {
    var matches = new List<int>();
    for (int i = 0; i < rows; i++) {
      if (condition(viewRow(i))) {
        matches.add(i);
      }
    }

    return viewSelection(matches, null);
  }

  /**
   * Returns a new selection view.
   */
  Matrix2d viewSelectionLike(List<int> rowOffsets, List<int> columnOffsets);

  /**
   * Sorts the matrix rows into ascending order, according to the natural ordering of the matrix
   * values in the given column.
   */
  Matrix2d viewSorted(int column) {
    throw 'Not implemented yet!';
    // return cern.colt.matrix.doublealgo.Sorting.mergeSort.sort(this,column);
  }

  /**
   * Constructs and returns a new stride view which is a sub matrix consisting of every i-th cell.
   */
  Matrix2d viewStrides(int rowStride, int columnStride) {
    return view()
        ..vStrides(rowStride, columnStride);
  }

  /**
   * Applies a procedure to each cell's value.  Iterates downwards from [rows - 1, columns - 1] to
   * [0, 0].
   */
  bool every(procedure) {
    for (int row = rows - 1; row >= 0; row--) {
      for (int column = columns - 1; column >= 0; column--) {
        if (!procedure(get(row, column))) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * 8 neighbor convolution. For efficient finite difference operations.  Applies a function to a
   * moving 3 x 3 window.  Does nothing if [rows] < 3 || [columns] < 3.
   */
  void convolve3x3(Matrix2d result, function) {
    if (function == null) {
      throw new ArgumentError('function must not be null.');
    }
    _checkValidShape(result);
    if (rows < 3 || columns < 3) {
      return;
    }

    int r = rows - 1;
    int c = columns - 1;
    double a00, a01, a02;
    double a10, a11, a12;
    double a20, a21, a22;
    for (int i = 1; i < r; i++) {
      a00 = get(i-1, 0);
      a01 = get(i-1, 1);

      a10 = get(i, 0);
      a11 = get(i, 1);

      a20 = get(i + 1, 0);
      a21 = get(i + 1, 1);

      for (int j = 1; j < c; j++) {
        // In each step six cells can be remembered in registers - they don't need to be reread from
        // slow memory in each step 3 instead of 9 cells need to be read from memory.
        a02 = get(i - 1, j + 1);
        a12 = get(i, j + 1);
        a22 = get(i + 1, j + 1);

        result.set(i,j, function(
          a00, a01, a02,
          a10, a11, a12,
          a20, a21, a22));

        a00 = a01;
        a10 = a11;
        a20 = a21;

        a01 = a02;
        a11 = a12;
        a21 = a22;
      }
    }
  }

  /**
   * Linear algebraic matrix-vector multiplication; result = this * y + result (if 'result' is
   * non-null).  Returns the result of the multiplication.
   */
  Matrix1d mul1d(Matrix1d y, Matrix1d result) {
    return mul1dWithParams(y, result, 1, (result == null ? 1 : 0), false);
  }

  /**
   * Linear algebraic matrix-vector multiplication; result = alpha * this * y + beta * result.
   */
  Matrix1d mul1dWithParams(Matrix1d y, Matrix1d result, double alpha, double beta,
      bool transposeA) {
    if (transposeA) {
      return viewTranspose().mul1dWithParams(y, result, alpha, beta, false);
    }
    if (result == null) {
      result = new DenseMatrix1d(_rows);
    }
    if (_columns != y.size || _rows > result.size) {
      throw new ArgumentError('Incompatible args: ${toStringShort()}, ${y.toStringShort()}, ' +
          '${result.toStringShort()}');
    }
    for (int i = rows - 1; i >= 0; i--) {
      double s = 0;
      for (int j = columns - 1; j >= 0; j--) {
        s += get(i,j) * y[j];
      }
      result[i] = alpha * s + beta * result[i];
    }
    return result;
  }

  /**
   * Linear algebraic matrix-matrix multiplication; result = this x b + b (if b is non-null).  The
   * result is returned.
   */
  Matrix2d mul2d(Matrix2d b, Matrix2d result) {
    return mul2dWithParams(b, result, 1.0, (result == null ? 1.0 : 0.0), false, false);
  }

  /**
   * Linear algebraic matrix-matrix multiplication; result = alpha * this x b + beta * result.
   */
  Matrix2d mul2dWithParams(Matrix2d b, Matrix2d result, double alpha, double beta,
      bool transposeThis, bool transposeB) {
    if (transposeThis) {
      return viewTranspose().mul2dWithParams(b, result, alpha, beta, false, transposeB);
    }
    if (transposeB) {
      return mul2dWithParams(b.viewTranspose(), result, alpha, beta, transposeThis, false);
    }

    int m = _rows;
    int n = _columns;
    int p = b._columns;

    if (result == null) {
      // TODO(klochek): This looks funny--shouldn't this be a 'like'?
      result = new DenseMatrix2d.withDimensions(m, p);
    }
    if (b._rows != n) {
      throw new ArgumentError('Matrix2D inner dimensions must agree:${toStringShort()}' +
          ' ${b.toStringShort()}');
    }
    if (result._rows != m || result._columns != p) {
      throw new ArgumentError('Incompatible result matrix: ${toStringShort()}, ' +
          '${b.toStringShort()}, ${result.toStringShort()}');
    }
    if (this == result || b == result) {
      throw new ArgumentError('Matrices must not be identical');
    }

    for (int j = p - 1; j >= 0; j--) {
      for (int i = m - 1; i >= 0; i--) {
        double s = 0.0;
        for (int k = n - 1; k >= 0; k--) {
          s += get(i, k) * b.get(k, j);
        }
        result.set(i, j, alpha * s + beta * result.get(i, j));
      }
    }
    return result;
  }

  /**
   * Returns the sum of all cells.
   */
  double sum() {
    if (size == 0) {
      return 0.0;
    }
    return aggregate((x) => x);
  }
}
