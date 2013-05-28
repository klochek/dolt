/*
 * Copyright � 1999 CERN - European Organization for Nuclear Research.
 * Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
 * is hereby granted without fee, provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear in supporting documentation.
 * CERN makes no representations about the suitability of this software for any purpose.
 * It is provided "as is" without expressed or implied warranty.
 *
 * Dart port Copyright (c) 2013 Chris Klochek.
 */

part of dolt;


/**
 * Dense 1-d matrix.  Internally holds one single contigous one-dimensional array.
 */
class DenseMatrix1d extends Matrix1d {
  /**
    * The elements of this matrix.
    */
  List<double> _elements;

  /**
   * Constructs a matrix with a copy of the given values.
   */
  DenseMatrix1d(List<double> values) {
    _setUp(values.length);
    _elements = new List<double>(size);
    assign(values);
  }

  /**
   * Constructs a matrix with [size] number of cells.
   */
  DenseMatrix1d.withSize(int size) {
    _setUp(size);
    _elements = new List<double>(size);
  }

  /**
   * Constructs a matrix view with the given parameters.
   */
  DenseMatrix1d.withParams(int size, List<double> elements, int zero, int stride) {
    _setUp(size, zero, stride);
    _elements = elements;
    _isView = true;
  }

  /**
   * Sets all cells to the state specified by value.
   */
  void assign(dynamic value) {
    if (value is List<double>) {
      if (!_isView) {
        if (value.length != size) {
          throw new ArgumentError('Lengths must be equal: length=${value.length} size=$size');
        }
        _elements.setAll(0, value);
      } else {
        super.assign(value);
      }
    } else if (value is double) {
      int index = this.index(0);
      int s = _stride;
      for (int i = size - 1; i >= 0; i--) {
        _elements[index] = value;
        index += s;
      }
    } else if (value is Function) {
      int s = _stride;
      int i = index(0);
      if (_elements == null) {
        throw 'Null elements error!';
      }

      // TODO(klochek): add function specialization optimization.
      for (int k = size - 1; k >= 0; k--) {
        _elements[i] = value(_elements[i]);
        i += s;
      }
    } else if (value is Matrix1d) {
      // Overriden for performance only.
      if (value is! DenseMatrix1d) {
        super.assign(value);
        return;
      }
      DenseMatrix1d other = value;
      if (other == this) {
        return;
      }
      _checkValidSize(other);
      if (!_isView && !other._isView) {
        _elements.setAll(0, other._elements);
      }
      if (haveSharedCells(other)) {
        Matrix1d c = other.copy();
        if (c is! DenseMatrix1d) {
          // Should not happen.
          return super.assign(value);
        }
        other = c;
      }

      List<double> elems = _elements;
      List<double> otherElems = other._elements;
      if (elems == null || otherElems == null) {
        throw 'Null elements error!';
      }
      int s = _stride;
      int ys = other._stride;

      int index = this.index(0);
      int otherIndex = other.index(0);
      for (int k = size - 1; k >= 0; k--) {
        elems[index] = otherElems[otherIndex];
        index += s;
        otherIndex += ys;
      }
    }
  }

  /**
   * Assigns the result of a function [f] to each cell in [this], so that this[i] =
   * f(this[i], y[i]).
   */
  assignWithFunc(Matrix1d y, function(double x, double y)) {
    // overriden for performance only
    if (y is! DenseMatrix1d) {
      super.assignWithFunc(y, function);
      return;
    }
    DenseMatrix1d other = y;
    _checkValidSize(y);
    List<double> elems = _elements;
    List<double> otherElems = other._elements;
    if (elems == null || otherElems == null) {
      throw 'Null elements error!';
    }
    int s = _stride;
    int ys = other._stride;

    int index = this.index(0);
    int otherIndex = other.index(0);

    // TODO(klochek): add function specialization optimizations.
    for (int k = size; k >= 0; k--) {
      elems[index] = function(elems[index], otherElems[otherIndex]);
      index += s;
      otherIndex += ys;
    }
  }

  /**
   * Returns the number of cells having non-zero values, up to an optional maximum of
   * [maxCardinality]; ignores tolerance.
   */
  int cardinality([int maxCardinality = -1]) {
    maxCardinality = maxCardinality > -1 ? maxCardinality : _size;
    int cardinality = 0;
    int index = this.index(0);
    int s = _stride;
    int i = _size - 1;
    while (i-- >= 0 && cardinality < maxCardinality) {
      if (_elements[index] != 0) {
        cardinality++;
      }
      index += s;
    }
    return cardinality;
  }

  /**
   * Returns the [double] value in the matrix cell value at coordinate [index].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  double operator[](int index) {
    return _elements[_zero + index * _stride];
  }

  /**
   * Returns [true] if both matrices share at least one identical cell.
   */
  bool haveSharedCellsRaw(Matrix1d other) {
    if (other is SelectedDenseMatrix1d) {
      return _elements == other._elements;
    } else if (other is DenseMatrix1d) {
      return _elements == other._elements;
    }
    return false;
  }

  /**
   * Returns the position of the element with the given relative [rank] within this matrix.
   */
  int index(int rank) {
    return _zero + rank * _stride;
  }

  /**
   * Construct and returns a new empty matrix of the same dynamic type as this matrix, with an
   * optionally-different specified [size].
   */
  Matrix1d like([int size = -1]) {
    size = size > -1 ? size : _size;
    return new DenseMatrix1d.withSize(size);
  }

  /**
   * Construct and returns a new 2-d matrix of the same dynamic type as this matrix, with specified
   * [rows] and [columns] sizes.
   * TODO(klochek): This is waiting on the Matrix2d implementation.
   */
  /*Matrix2d like2D(int rows, int columns) {
    return new DenseMatrix2d(rows, columns);
  }*/

  /**
   * Sets the [double] value in the matrix cell value at coordinate [index].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  void operator[]=(int index, double value) {
    _elements[_zero + index * _stride] = value;
  }

  /**
   * Swaps each element this[i] with other[i].
   */
  void swap(Matrix1d other) {
    // overriden for performance only
    if (other is! DenseMatrix1d) {
      super.swap(other);
    }
    DenseMatrix1d y = other as DenseMatrix1d;
    if (y == this) {
      return;
    }
    _checkValidSize(y);

    List<double> elems = _elements;
    List<double> otherElems = y._elements;
    if (elems == null || otherElems == null) {
      throw 'Null elements';
    }
    int s = _stride;
    int ys = y._stride;
    int index = this.index(0);

    int otherIndex = y.index(0);
    for (int k = size - 1; k >= 0; k--) {
      double tmp = elems[index];
      elems[index] = otherElems[otherIndex];
      otherElems[otherIndex] = tmp;
      index += s;
      otherIndex += ys;
    }
  }

  /**
   * Fills the cell values in the optionally specified list [values], returning it.  If a list is
   * specified, it must be of the same size as this matrix.
   */
  void toArray([List<double> values]) {
    if (values != null) {
      if (values.length < size) {
        throw new ArgumentError("values too small");
      }
    } else {
      values = new List<double>(size);
    }
    if (!_isView) {
      values.setAll(0, _elements);
    } else {
      super.toArray(values);
    }
  }

  /**
   * Construct and returns a new selection view.
   */
  Matrix1d viewSelectionLike(List<int> offsets) {
    return new SelectedDenseMatrix1d(_elements, offsets);
  }

  /**
   * Returns the dot product this matrix and [y];
   */
  double dotRange(Matrix1d y, int from, int length) {
    if (y is! DenseMatrix1d) {
      return super.dot(y, from, length);
    }
    DenseMatrix1d yy = y;

    int tail = from + length;
    if (from < 0 || length < 0) {
      return 0.0;
    }
    if (size < tail) {
      tail = size;
    }
    if (y.size < tail) {
      tail = y.size;
    }
    int min = tail - from;

    int i = index(from);
    int j = yy.index(from);
    int s = _stride;
    int ys = yy._stride;
    List<double> elems = _elements;
    List<double> yElems = yy._elements;
    if (elems == null || yElems == null) {
      throw 'Null elements!';
    }

    double sum = 0.0;
    // Loop unrolling.
    i -= s;
    j -= ys;
    for (int k = (min / 4).floor() - 1; k >= 0; k--) {
      sum += elems[i += s] * yElems[j += ys] +
        elems[i += s] * yElems[j += ys] +
        elems[i += s] * yElems[j += ys] +
        elems[i += s] * yElems[j += ys];
    }
    for (int k = (min % 4) - 1; k >= 0; k--) {
      sum += elems[i += s] * yElems[j += ys];
    }
    return sum;
  }

  /**
   * Returns the sum of all cells.
   */
  double sum() {
    double sum = 0.0;
    int s=stride;
    int i=index(0);
    List<double> elems = _elements;
    if (elems == null) {
      throw 'Null elements!';
    }
    for (int k = size - 1; k >= 0; k--) {
      sum += elems[i];
      i += s;
    }
    return sum;
  }
}
