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
 * Sparse hashed 2-d matrix holding double elements. Cells that are never set to non-zero values do
 * not use any memory.  Cells that switch from zero to non-zero state do use memory.  Cells that
 * switch back from non-zero to zero state also do use memory.  This class offers constant
 * expected-time complexity for the basic operations.
 */
class SparseMatrix2d extends Matrix2d {
  /// The elements of the matrix.
  Map<int, double> elements;

  /**
   * Constructs a matrix with a copy of the given values.  [values] is required to have the form
   * values[row][column], and have exactly the same number of columns in every row.
   * The values are copied, so subsequent changes in [values] are not reflected in the matrix, and
   * vice-versa.
   */
  SparseMatrix2d(List<List<double>> values) {
    _setUp(values.length, values.length == 0 ? 0 : values[0].length);
    elements = new Map<int, double>();
    assign(values);
  }

  /**
   * Constructs a matrix with a given number of rows and columns.  All entries are initially 0.
   */
  SparseMatrix2d.withDimensions(int rows, int columns) {
    elements = new Map<int, double>();
    _setUp(rows, columns);
  }

  /**
   * Constructs a view with the given parameters.
   */
  SparseMatrix2d.withValues(int rows, int columns, Map<int, double> elements, int rowZero,
      int columnZero, int rowStride, int columnStride) {
    _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);
    this.elements = elements;
    _isView = true;
  }

  /**
   * Sets all cells to the state specified by [value].
   */
  void assign(dynamic value) {
    // Overriden for performance only.
    if (value is Function) {
      // TODO(klochek): specialize for performance, later.
      /*if (!isView && function instanceof cern.jet.math.Mult) {
        this.elements.assign(function);
      }*/
      super.assign(function);
    } else if (value is Matrix2d) {
      if (source is! SparseMatrix2d) {
        super.assign(source);
        return;
      }
      SparseMatrix2d other = source;
      if (other == this) {
        return;
      }
      checkShape(other);

      if (!isView && !other.isView) {
        elements.assign(other.elements);
        return;
      }
    } else {
      if (!isView && value == 0) {
        elements.clear();
      }
    }
    super.assign(value);
  }

  void assignWithFunc(Matrix2d y, function) {
    if (isView) {
      super.assign(y, function);
      return;
    }

    checkShape(y);

    // TODO(klochek): Specialize for performance, later.
    super.assign(y, function);
  }

  /**
   * Returns the number of cells having non-zero values.
   */
  int cardinality() {
    if (!isView) {
      return elements.length;
    } else {
      return super.cardinality();
    }
  }

  void ensureCapacity(int minCapacity) {
    // TODO(klochek): Probably useless.
  }

  Matrix2d forEachNonZero(function) {
    if (!isView) {
      elements.forEach((k, v) {
        int i = k / columns;
        int j = k % columns;
        double r = function(i, j, v);
        if (r != v) {
          elements[k] = v;
        }
      });
    } else {
      super.forEachNonZero(function);
    }
  }

  /**
   * Returns the matrix cell value at coordinate this[row,column].
   */
  double get(int row, int column) {
    return elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride];
  }

  /**
   * Returns [true] if both matrices share common cells.  More formally, returns [true] if at least
   * one of the following conditions is met: the receiver is a view of the other matrix; the other
   * matrix is a view of the receiver; this == other.
   */
  bool haveSharedCellsRaw(Matrix2d other) {
    if (other is SelectedSparseMatrix2d) {
      return elements == other.elements;
    } else if (other is SparseMatrix2d) {
      return elements == other.elements;
    }
    return false;
  }

  /**
   * Returns the position of the given coordinate within the (virtual or non-virtual) internal
   * 1-dimensional array.
   */
  int index(int row, int column) {
    return _rowZero + row * _rowStride + _columnZero + column * _columnStride;
  }

  /**
   * Construct and returns a new empty matrix of the same dynamic type as the receiver, having
   * the specified number of rows and columns.
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
   * Construct and returns a new 1-d matrix of the corresponding dynamic type, sharing the
   * same cells.
   */
  Matrix1d like1dWithDimensions(int size, int offset, int stride) {
    return new SparseMatrix1d.withOffsets(size, elements, offset, stride);
  }

  /**
   * Sets the matrix cell at coordinate [row,column] to the specified value.
   */
  void set(int row, int column, double value) {
    int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;

    if (value == 0) {
      elements.remove(index);
    } else {
      elements[index] = value;
    }
  }

  void trimToSize() {
    // TODO(klochek): Useless?
  }

  /**
   * Construct and returns a new selection view.
   */
  Matrix2d viewSelectionLike(List<int> rowOffsets, List<int> columnOffsets) {
    return new SelectedSparseMatrix2d(elements, rowOffsets, columnOffsets, 0);
  }

  Matrix1d mul1dWithParams(Matrix1d y, Matrix1d z, double alpha, double beta, bool transposeA) {
    int m = rows;
    int n = columns;
    if (transposeA) {
      m = columns;
      n = rows;
    }

    bool ignore = (z == null);
    if (z == null) {
      z = new DenseMatrix1d(m);
    }

    if (!(!isView && y is DenseMatrix1d && z is DenseMatrix1d)) {
      return super.mul2dWithParams(y, z, alpha, beta, transposeA);
    }

    if (n != y.size() || m > z.size()) {
      throw new ArgumentError('Incompatible args: ' +
          ((transposeA ? viewDice() : this).toStringShort()) + ', ' + y.toStringShort() + ', ' +
          z.toStringShort());
    }

    if (!ignore) {
      var c = beta / alpha;
      z.assign((v) => v * c);
    }

    DenseMatrix1d zz = z;
    List<double> zElements = zz.elements;
    int zStride = zz.stride;
    int zi = z.index(0);

    DenseMatrix1d yy = y;
    List<double> yElements = yy.elements;
    int yStride = yy.stride;
    int yi = y.index(0);

    if (yElements == null || zElements == null) {
      throw 'Null elements should never happen';
    }

    elements.forEach((k, v) {
      int i = k / columns;
      int j = k % columns;
      if (transposeA) {
        int tmp = i;
        i = j;
        j = tmp;
      }
      zElements[zi + zStride * i] += value * yElements[yi + yStride * j];
    });

    if (alpha != 1) {
      z.assign((v) => v * alpha);
    }
    return z;
  }

  Matrix2d mul2dWithParams(Matrix2d B, Matrix2d C, double alpha, double beta, bool transposeA,
      bool transposeB) {

    if (isView) {
      return super.mul2dWithParams(B, C, alpha, beta, transposeA, transposeB);
    }

    if (transposeB) {
      B = B.viewTranspose();
    }
    int m = rows;
    int n = columns;
    if (transposeA) {
      m = columns;
      n = rows;
    }
    int p = B.columns;
    bool ignore = (C == null);
    if (C == null) {
      C = new DenseMatrix2d.withDimensions(m, p);
    }

    if (B.rows != n) {
      throw new ArgumentError('Matrix2D inner dimensions must agree:' + toStringShort() + ', ' +
          (transposeB ? B.viewTranspose() : B).toStringShort());
    }
    if (C.rows != m || C.columns != p) {
      throw new ArgumentError('Incompatibel result matrix: ' + toStringShort() + ', ' +
          (transposeB ? B.viewTranspose() : B).toStringShort() + ', ' + C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError('Matrices must not be identical');
    }

    if (!ignore) {
      C.assign((v) => v * beta);
    }

    // cache views
    List<Matrix1d> Brows = new List<Matrix1d>(n);
    for (int i = n - 1; i >= 0; i--) {
      Brows[i] = B.viewRow(i);
    }

    List<Matrix1d> Crows = new List<Matrix1d>(m);
    for (int i = m - 1; i >= 0; i--) {
      Crows[i] = C.viewRow(i);
    }

    elements.forEach((k, v) {
      int i = (k / columns).floor();
      int j = k % columns;
      var multiplicator = v * alpha;
      var func = (a, b) => a + b * multiplicator;
      if (!transposeA) {
        Crows[i].assignWithFunc(Brows[j], func);
      } else {
        Crows[j].assignWithFunc(Brows[i], func);
      }
    });

    return C;
  }
}
