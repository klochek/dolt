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
 * Selection view on sparse 2-d matrices holding double elements.
 */
class SelectedSparseMatrix2d extends Matrix2d {
  /// The elements of the matrix.
  Map<int, double> elements;

  /// The row offsets of the visible cells of this matrix.
  List<int> _rowOffsets;

  /// The column offsets of the visible cells of this matrix.
  List<int> _columnOffsets;

  int _offset;

  /**
   * Constructs a matrix view with the given parameters.
   */
  SelectedSparseMatrix2d.withParams(int rows, int columns, Map<int, double> elements, int rowZero,
      int columnZero, int rowStride, int columnStride, List<int> rowOffsets,
      List<int> columnOffsets, int offset) {
    _setUpWithParams(rows, columns, rowZero, columnZero, rowStride, columnStride);

    this.elements = elements;
    _rowOffsets = rowOffsets;
    _columnOffsets = columnOffsets;
    _offset = offset;

    _isView = true;
  }

  /**
   * Constructs a matrix view with the given parameters.
   */
  SelectedSparseMatrix2d(Map<int, double> elements, List<int> rowOffsets, List<int> columnOffsets,
      int offset) :
        this.withParams(rowOffsets.length, columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets,
            columnOffsets, offset);

  /**
   * Returns the position of the given absolute rank within the (virtual or non-virtual) internal
   * 1-dimensional array.
   */
  int _columnOffset(int absRank) {
    return _columnOffsets[absRank];
  }

  /**
   * Returns the position of the given absolute rank within the (virtual or non-virtual) internal
   * 1-dimensional array.
   */
  int _rowOffset(int absRank) {
    return _rowOffsets[absRank];
  }

  /**
   * Returns the matrix cell value at coordinate [row, column].
   */
  double get(int row, int column) {
    return elements[_offset + _rowOffsets[_rowZero + row * _rowStride] +
        _columnOffsets[_columnZero + column * _columnStride]];
  }

  /**
   * Returns [true] if both matrices share common cells.  More formally, returns [true] if at least
   * one of the following conditions is met: the receiver is a view of the other matrix; the other
   * matrix is a view of the receiver; this == other.
   */
  bool haveSharedCellsRaw(Matrix2d other) {
    if (other is SelectedSparseMatrix2d || other is SparseMatrix2d) {
      return elements == other.elements;
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
   * Construct and returns a new empty matrix of the same dynamic type as the receiver, having the
   * specified number of rows and columns.
   */
  Matrix2d likeWithDimensions(int rows, int columns) {
    return new SparseMatrix2d(rows, columns);
  }

  /**
   * Construct and returns a new 1-d matrix of the corresponding dynamic type, entirelly independent
   * of the receiver.
   */
  Matrix1d like1d(int size) {
    return new SparseMatrix1d(size);
  }

  /**
   * Construct and returns a new 1-d matrix of the corresponding dynamic type, sharing the same
   * cells.
   */
  Matrix1d like1dWithParams(int size, int zero, int stride) {
    throw 'This should not happen';
  }

  /**
   * Sets the matrix cell at coordinate [row,column] to the specified value.
   */
  void set(int row, int column, double value) {
    int index = _offset + _rowOffsets[_rowZero + row * _rowStride] +
        _columnOffsets[_columnZero + column * _columnStride];

    if (value == 0) {
      elements.remove(index);
    } else {
      elements[index] = value;
    }
  }

  /**
   * Sets up a matrix with a given number of rows and columns.
   */
  void _setUp(int rows, int columns) {
    super._setUp(rows, columns);
    _rowStride = 1;
    _columnStride = 1;
    _offset = 0;
  }

  /**
   * Self modifying version of [viewTranspose].
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
   * ([viewPart(...)]), then apply this method to the sub-range view.
   */
  Matrix1d viewColumn(int column) {
    checkColumn(column);
    int viewSize = _rows;
    int viewZero = _rowZero;
    int viewStride = _rowStride;
    List<int> viewOffsets = _rowOffsets;
    int viewOffset = _offset + _columnOffset(_columnRank(column));
    return new SelectedSparseMatrix1d(viewSize, elements, viewZero, viewStride, viewOffsets,
        viewOffset);
  }

  /**
   * Constructs and returns a new slice view representing the columns of the given row.  The
   * returned view is backed by this matrix, so changes in the returned view are reflected in this
   * matrix, and vice-versa.  To obtain a slice view on subranges, construct a sub-ranging view
   * ([viewPart(...)]), then apply this method to the sub-range view.
   */
  Matrix1d viewRow(int row) {
    checkRow(row);
    int viewSize = _columns;
    int viewZero = _columnZero;
    int viewStride = _columnStride;
    List<int> viewOffsets = _columnOffsets;
    int viewOffset = _offset + _rowOffset(_rowRank(row));
    return new SelectedSparseMatrix1d(viewSize, elements, viewZero, viewStride, viewOffsets,
        viewOffset);
  }

  /**
   * Construct and returns a new selection view.
   */
  Matrix2d viewSelectionLike(List<int> rowOffsets, List<int> columnOffsets) {
    return new SelectedSparseMatrix2d(elements, rowOffsets, columnOffsets, _offset);
  }
}
