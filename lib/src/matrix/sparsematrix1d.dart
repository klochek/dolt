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
 * Sparse hashed 1-d matrix.  Cells that are never set to non-zero values do not use any memory.
 * Cells that switch from zero to non-zero state do use memory.  Cells that switch back from
 * non-zero to zero state will be garbage collected.
 * This class offers constant expected time complexity for the basic get/set operations.
 */
class SparseMatrix1d extends Matrix1d {
  /// The elements of the matrix.
  Map<int, double> _elements;

  /**
   * Constructs a matrix with a copy of the given [values].
   */
  SparseMatrix1d(List<double> values) :
      this.withSize(values.length, values);

  /**
   * Constructs a matrix with a given [size].  All entries are initially 0.
   */
  SparseMatrix1d.withSize(int size, [List<double> values]) :
      this.withParams(size, (size / 1000).floor(), 0.2, 0.5, values);

  /**
   * Constructs a matrix with a given number of parameters.  All entries are initially 0.
   */
  SparseMatrix1d.withParams(int size, int initialCapacity, double minLoadFactor,
      double maxLoadFactor, [List<double> values]) {
    _setUp(size);
    _elements = new Map<int, double>();
    assign(values);
    // TODO(klochek): Investigate efficiency of their hashmap, and possibly deprecate this ctor.
    //_elements = new OpenIntDoubleHashMap(initialCapacity, minLoadFactor, maxLoadFactor);
  }

  /**
   * Constructs a matrix view with a given number of parameters.
   */
  SparseMatrix1d.withOffsets(int size, Map<int, double> elements, int offset, int stride) {
    _setUp(size, offset, stride);
    _elements = elements;
    _isView = true;
  }

  /**
   * Sets all cells to the state specified by [value].
   */
  void assign(dynamic value) {
    // overriden for performance only
    if (value is double && !_isView && value == 0) {
      _elements.clear();
    } else {
      super.assign(value);
    }
  }

  /**
   * Returns the number of cells having non-zero values.
   */
  int cardinality([maxCardinality = -1]) {
    if (!_isView) {
      return _elements.length;
    } else {
      return super.cardinality(maxCardinality);
    }
  }

  /**
   * Ensures that the receiver can hold at least the specified number of non-zero cells without
   * needing to allocate new internal memory.
   */
  void ensureCapacity(int minCapacity) {
    // TODO(klochek): maybe deprecate?
    throw "Not implemented yet!";
    //this.elements.ensureCapacity(minCapacity);
  }

  /**
   * Returns the matrix cell value at coordinate [index].
   */
  double operator[](int index) {
    var el = _elements[_zero + index * _stride];
    return el == null ? 0.0 : el;
  }

  /**
   * Returns [true] if both this matrix and [other] share at least one identical cell.
   */
  bool haveSharedCellsRaw(Matrix1d other) {
    // TODO(klochek): implement SelectedSparseMatrix1d.
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
    // overriden for manual inlining only
    return _zero + rank * _stride;
  }

  /**
   * Construct and returns a new empty matrix of the same dynamic type as this matrix, with the
   * specified [size].
   */
  Matrix1d like([int size = -1]) {
    size = size > 0 ? size : _size;
    return new SparseMatrix1d(size);
  }

  /**
   * Construct and returns a new sparse 2-d matrix with specified [rows] and [columns].
   */
  /*Matrix2d like2D(int rows, int columns) {
    return new SparseMatrix2d(rows,columns);
  }*/

  /**
   * Sets the [double] value in the matrix cell value at coordinate [index].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  operator[]=(int index, double value) {
    int i = _zero + index * _stride;
    if (value == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = value;
    }
  }

  /**
   * Releases any superfluous memory created by explicitly putting zero values into cells formerly
   * having non-zero values.
   */
  void trimToSize() {
    // TODO(klochek): might be pointless, but worth investigating.
    // this.elements.trimToSize();
  }

  /**
   * Construct and returns a new selection view, given the provided [offsets].
   */
  Matrix1d viewSelectionLike(List<int> offsets) {
    return new SelectedSparseMatrix1d(_elements, offsets);
  }
}
