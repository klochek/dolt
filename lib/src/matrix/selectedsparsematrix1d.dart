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
 * Selection view on sparse 1-d matrices.
 */
class SelectedSparseMatrix1d extends Matrix1d {

  /// The elements of the matrix.
  Map<int, double> _elements;

  /// The offsets of visible indexes of this matrix.
  List<int> _offsets;

  /// The offset.
  int _actualOffset;

  /**
   * Constructs a matrix view with the given parameters.
   */
  SelectedSparseMatrix1d(int size, Map<int, double> elements, int zero, int stride,
      List<int> offsets, int offset) {
    _setUp(size, zero, stride);

    _elements = elements;
    _offsets = offsets;
    _actualOffset = offset;
    _isView = true;
  }

  /**
   * Constructs a matrix view with the given parameters.
   */
  SelectedSparseMatrix1d.fromElements(Map<int, double> elements, List<int> offsets) :
    this(offsets.length, elements, 0, 1, offsets, 0);

  /**
   * Returns the [int] position of the given absolute rank within the (virtual or non-virtual)
   * internal 1-dimensional array.
   */
  int _offset(int absRank) {
    return _offsets[absRank];
  }

  /**
   * Returns the [double] value in the matrix cell value at coordinate [index].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  double operator[](int index) {
    return _elements[_actualOffset + _offsets[_zero + index * _stride]];
  }

  /**
   * Sets the [double] value in the matrix cell value at coordinate [index].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  void operator[]=(int index, double value) {
    int i = _actualOffset + _offsets[_zero + index * _stride];
    if (value == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = value;
    }
  }

  /**
   * Returns [true] if both this matrix and [other] share at least one identical cell.
   */
  bool haveSharedCellsRaw(Matrix1d other) {
    if (other is SelectedSparseMatrix1d) {
      return _elements == other._elements;
    } else if (other is SparseMatrix1d) {
      return _elements == other._elements;
    }
    return false;
  }

  /**
   * Returns the position of the element with the given relative [rank] within this matrix.
   */
  int index(int rank) {
    return _actualOffset + _offsets[_zero + rank * _stride];
  }

  /**
   * Construct and returns a new empty matrix of the same dynamic type as this matrix, with the
   * specified [size].
   */
  Matrix1d like([int size = -1]) {
    size = size > -1 ? size : _size;
    return new SparseMatrix1d(size);
  }

  /**
   * Construct and returns a new sparse 2-d matrix with specified [rows] and [columns].
   */
  /*Matrix2d like2D(int rows, int columns) {
    return new SparseMatrix2d(rows, columns);
  }*/

  /**
   * Sets up a matrix with a given [size].
   */
  void _setUp(int size, [int zero = 0, int stride = 1]) {
    // We ignore 'zero; and 'stride' for selection matrices.
    super._setUp(size);
    _stride = 1;
    _actualOffset = 0;
  }

  /**
   * Construct and returns a new selection view with the supplied [offsets].
   */
  Matrix1d viewSelectionLike(List<int> offsets) {
    return new SelectedSparseMatrix1d(_elements, _offsets);
  }
}
