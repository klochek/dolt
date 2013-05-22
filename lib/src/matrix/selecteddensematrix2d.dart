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
 * Selection view on dense 2-d matrices.
 */
class SelectedDenseMatrix2d extends Matrix2d {
  /// The elements of this matrix.
  List<double> _elements;

  /// The row offsets of the visible cells of this matrix.
  List<int> _rowOffsets;

  /// The column offsets of the visible cells of this matrix.
  List<int> _columnOffsets;

  /// The offset.
  int _offset;

  /**
   * Constructs a matrix view with the given parameters.
   */
  SelectedDenseMatrix2d(List<double> elements, List<int> rowOffsets, List<int> columnOffsets,
      int offset)
      : this._view(rowOffsets.length, columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets,
          columnOffsets, offset);

  /**
   * Constructs a matrix view with the given parameters.
   */
  SelectedDenseMatrix2d._view(int rows, int columns, List<double> elements, int rowZero,
      int columnZero, int rowStride, int columnStride, List<int> rowOffsets,
      List<int> columnOffsets, int offset) {
    _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);

    _elements = elements;
    _rowOffsets = rowOffsets;
    _columnOffsets = columnOffsets;
    _offset = offset;

    _isView = true;
  }

  /**
   * Returns the [int] position of the given absolute rank within the (virtual or non-virtual)
   * internal 1-dimensional array.
   */
  int _columnOffset(int absRank) {
    return _columnOffsets[absRank];
  }

  /**
   * Returns the [int] position of the given absolute rank within the (virtual or non-virtual)
   * internal 1-dimensional array.
   */
  int _rowOffset(int absRank) {
    return _rowOffsets[absRank];
  }

  /**
   * Returns the [double] value in the matrix cell value at [row] and [column].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  double get(int row, int column) {
    return _elements[_offset + _rowOffsets[_rowZero + row * _rowStride] +
        _columnOffsets[_columnZero + column * _columnStride]];
  }

  /**
   * Returns [true] if both matrices share at least one identical cell.
   */
  bool haveSharedCellsRaw(Matrix2d other) {
    if (other is SelectedDenseMatrix2d) {
      return _elements == other._elements;
    } else if (other is DenseMatrix2d) {
      return _elements == other._elements;
    }
    return false;
  }

  /**
   * Returns the position of the given coordinate within the (virtual or non-virtual) internal
   * 1-dimensional array.
   */
  int index(int row, int column) {
    return _offset + _rowOffsets[_rowZero + row * _rowStride] +
        _columnOffsets[_columnZero + column * _columnStride];
  }

  /**
   * Construct and returns a new empty matrix of the same dynamic type as this matrix, having the
   * same number of rows and columns.
   */
  Matrix2d like(int rows, int columns) {
    return new DenseMatrix2d(rows, columns);
  }

  /**
   * Construct and returns a new 1-d matrix of the corresponding dynamic type of this matrix, with
   * the specified [size].
   */
  Matrix1d like1d(int size) {
    return new DenseMatrix1d(size);
  }

  /**
   * Construct and returns a new 1-d matrix of the corresponding dynamic type of this matrix, with
   * specified [size], [zero], and [stride].
   */
  Matrix1d like1dWithDimensions(int size, int zero, int stride) {
    // TODO(klochek): This used to be protected; refactor so that a client cannot call this.
    throw 'Internal error'; // this method is never called since viewRow() and viewColumn are overridden properly.
  }

  /**
   * Sets the matrix cell at coordinate [row] [column] to the specified [value].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  void set(int row, int column, double value) {
    _elements[_offset + _rowOffsets[_rowZero + row * _rowStride] +
        _columnOffsets[_columnZero + column * _columnStride]] = value;
  }

  /**
   * Sets up a matrix with a given number of [rows] and [columns].
   */
  void _setUp(int rows, int columns) {
    super._setUp(rows, columns);
    _rowStride = 1;
    _columnStride = 1;
    _offset = 0;
  }

  /**
   * Self modifying version of viewTranspose().
   */
  void vTranspose() {
    super.vTranspose();
    List<int> tmp = _rowOffsets;
    _rowOffsets = _columnOffsets;
    _columnOffsets = tmp;
    _isView = true;
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
    int viewZero = _rowZero;
    int viewStride = _rowStride;
    List<int> viewOffsets = _rowOffsets;
    int viewOffset = _offset + _columnOffset(_columnRank(column));
    return new SelectedDenseMatrix1d(viewSize, _elements, viewZero, viewStride, viewOffsets,
        viewOffset);
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
    int viewZero = _columnZero;
    int viewStride = _columnStride;
    List<int> viewOffsets = _columnOffsets;
    int viewOffset = _offset + _rowOffset(_rowRank(row));
    return new SelectedDenseMatrix1d(viewSize, _elements, viewZero, viewStride, viewOffsets,
        viewOffset);
  }

  /**
   * Returns a new selection view.
   */
  Matrix2d viewSelectionLike(List<int> rowOffsets, List<int> columnOffsets) {
    return new SelectedDenseMatrix2d(_elements, _rowOffsets, _columnOffsets, _offset);
  }
}
