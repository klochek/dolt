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
 * Dense 2-d matrix.
 */
class DenseMatrix2d extends Matrix2d {

  /// Block-size to use when doing multiplication.
  int blockSize = 30000;

  /// The elements of this matrix.  Elements are stored in row major form
  /// (index== row * columns + column).
  List<double> _elements;

  /**
   * Constructs a matrix with a copy of the given values.
   */
  DenseMatrix2d(List<List<double>> values) {
    _setUp(values.length, values.length == 0 ? 0 : values[0].length);
    _isView = false;
    _elements = new List<double>.filled(_rows * _columns, 0.0);
    assign(values);
  }

  /**
   * Constructs a matrix with a given number of [rows] and [columns].  All entries are initialized
   * to 0.
   */
  DenseMatrix2d.withDimensions(int rows, int columns) {
    _setUp(rows, columns);
    _elements = new List<double>.filled(rows * columns, 0.0);
  }

  /**
   * Constructs a view with the given parameters.
   */
  DenseMatrix2d.withParams(int rows, int columns, List<double> elements, int rowZero,
      int columnZero, int rowStride, int columnStride) {
    _setUpWithParams(rows, columns, rowZero, columnZero, rowStride, columnStride);
    _elements = elements;
    _isView = true;
  }

  /**
   * Sets all cells to the state specified by [value].
   */
  void assign(dynamic value) {
    if (value is List<List<double>>) {
      List<List<double>> values = value;
      if (!_isView) {
        if (values.length != _rows) {
          throw new ArgumentError('Must have same number of rows: rows=${values.length}' +
              'rows()=$_rows');
        }
        int i = _columns * (_rows - 1);
        for (int row = rows - 1; row >= 0; row--) {
          var currentRow = values[row];
          if (currentRow.length != columns) {
            throw new ArgumentError('Must have same number of columns in every row: columns=' +
                '${currentRow.length} columns()=$_columns');
          }
          _elements.setRange(i, i + _columns, currentRow);
          i -= _columns;
        }
      } else {
        super.assign(values);
      }
    } else if (value is double) {
      int idx = index(0, 0);
      int cs = _columnStride;
      int rs = _rowStride;
      for (int row = rows - 1; row >= 0; row--) {
        for (int i = idx, column = columns - 1; column >= 0; column--) {
          _elements[i] = value;
          i += cs;
        }
        idx += rs;
      }
    } else if (value is Function) {
      if (_elements == null) {
        throw 'Elements must not be null.';
      }
      int idx = index(0, 0);
      int cs = _columnStride;
      int rs = _rowStride;

      // TODO(klochek): add specializations for performance.
      for (int row = rows - 1; row >= 0; row--) {
        for (int i = idx, column = columns - 1; column >= 0; column--) {
          _elements[i] = value(_elements[i]);
          i += cs;
        }
        idx += rs;
      }
    } else if (value is Matrix2d) {
      Matrix2d source = value;

      if (source is! DenseMatrix2d) {
        return super.assign(source);
      }
      DenseMatrix2d other = source as DenseMatrix2d;
      if (other == this) {
        return;
      }
      _checkValidShape(other);

      if (!_isView && !other._isView) {
        _elements.setAll(0, other._elements);
        return;
      }

      if (haveSharedCells(other)) {
        Matrix2d c = other.copy();
        if (c is! DenseMatrix2d) {
          return super.assign(other);
        }
        other = c;
      }

      var elems = _elements;
      var otherElems = other._elements;
      if (elems == null || otherElems == null) {
        throw 'Null elements should be impossible';
      }
      int cs = _columnStride;
      int ocs = other._columnStride;
      int rs = _rowStride;
      int ors = other._rowStride;

      int otherIndex = other.index(0, 0);
      int idx = index(0, 0);
      for (int row = rows - 1; row >= 0; row--) {
        for (int i = idx, j = otherIndex, column = columns - 1; column >= 0; column--) {
          elems[i] = otherElems[j];
          i += cs;
          j += ocs;
        }
        idx += rs;
        otherIndex += ors;
      }
    }
  }


  /**
   * Assigns the result of a function to each cell.
   */
  void assignWithFunc(Matrix2d y, function) {
    if (y is! DenseMatrix2d) {
      return super.assign(y, function);
    }
    DenseMatrix2d other = y as DenseMatrix2d;
    _checkValidShape(y);

    var elems = _elements;
    var otherElems = other._elements;
    if (elems == null || otherElems == null) {
      throw 'Null elements should be impossble';
    }
    int cs = _columnStride;
    int ocs = other._columnStride;
    int rs = _rowStride;
    int ors = other._rowStride;

    int otherIndex = other.index(0, 0);
    int idx = index(0, 0);

    // TODO(klochek): add specializations for performance.
    for (int row = rows - 1; row >= 0; row--) {
      for (int i = idx, j = otherIndex, column = columns - 1; column >= 0; column--) {
        elems[i] = function(elems[i], otherElems[j]);
        i += cs;
        j += ocs;
      }
      idx += rs;
      otherIndex += ors;
    }
  }

  /**
   * Returns the [double] value in the matrix cell value at [row] and [column].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  double get(int row, int column) {
    return _elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride];
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
    return _rowZero + row * _rowStride + _columnZero + column * _columnStride;
  }

  /**
   * Construct and returns a new empty matrix of the same dynamic type as this matrix, having the
   * specified number of [rows] and [columns].
   */
  Matrix2d likeWithDimensions(int rows, int columns) {
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
    return new DenseMatrix1d.withParams(size, _elements, zero, stride);
  }

  /**
   * Sets the matrix cell at coordinate [row] [column] to the specified [value].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  void set(int row, int column, double value) {
    _elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride] = value;
  }

  /**
   * Returns a new selection view.
   */
  Matrix2d viewSelectionLike(List<int> rowOffsets, List<int> columnOffsets) {
    return new SelectedDenseMatrix2d(_elements, rowOffsets, columnOffsets, 0);
  }

  /**
   * 8 neighbor convolution. For efficient finite difference operations.  Applies a function to a
   * moving 3 x 3 window.  Does nothing if [rows] < 3 || [columns] < 3.
   */
  void convolve3x3(Matrix2d b, function) {
    if (b is! DenseMatrix2d) {
      super.convolve3x3(b, function);
      return;
    }
    if (function == null) {
      throw new ArgumentError('function must not be null.');
    }
    _checkValidShape(b);
    int r = rows - 1;
    int c = columns - 1;
    if (rows < 3 || columns < 3) {
      return;
    }

    DenseMatrix2d bb = b as DenseMatrix2d;
    int this_rs = _rowStride;
    int b_rs = bb._rowStride;
    int this_cs = _columnStride;
    int b_cs = bb._columnStride;
    var elems = _elements;
    var b_elems = bb._elements;
    if (elems == null || b_elems == null) {
      throw 'Null elements should be impossible.';
    }

    int this_index = index(1, 1);
    int b_index = bb.index(1, 1);
    for (int i = 1; i < r; i++) {
      double this00, this01, this02;
      double this10, this11, this12;
      double this20, this21, this22;

      int b11 = b_index;

      int _this02 = this_index - this_rs - this_cs;
      int _this12 = _this02 + this_rs;
      int _this22 = _this12 + this_rs;

      this00 = elems[_this02]; _this02 += this_cs;   this01 = elems[_this02];
      this10 = elems[_this12]; _this12 += this_cs;   this11 = elems[_this12];
      this20 = elems[_this22]; _this22 += this_cs;   this21 = elems[_this22];

      for (int j = 1; j < c; j++) {
        // In each step 3 instead of 9 cells need to be read from memory.
        this02 = elems[_this02 += this_cs];
        this12 = elems[_this12 += this_cs];
        this22 = elems[_this22 += this_cs];

        b_elems[b11] = function(
          this00, this01, this02,
          this10, this11, this12,
          this20, this21, this22);
        b11 += b_cs;

        // Move remembered cells.
        this00 = this01; this01 = this02;
        this10 = this11; this11 = this12;
        this20 = this21; this21 = this22;
      }
      this_index += this_rs;
      b_index += this_rs;
    }
  }

  /**
   * Linear algebraic matrix-vector multiplication; result = alpha * this * y + beta * result.
   */
  Matrix1d mul1dWithParams(Matrix1d y, Matrix1d result, double alpha, double beta,
      bool transposeThis) {
    if (transposeThis) {
      return viewTranspose().mul1dWithParams(y, result, alpha, beta, false);
    }
    if (result == null) {
      result = new DenseMatrix1d(_rows);
    }
    if (!(y is DenseMatrix1d && result is DenseMatrix1d)) {
      return super.mul1dWithParams(y, result, alpha, beta, transposeThis);
    }

    if (_columns != y.size || _rows > result.size) {
      throw new ArgumentError('Incompatible args: ${toStringShort()}, ${y.toStringShort()}, ' +
          result.toStringShort());
    }
    DenseMatrix1d yy = y;
    DenseMatrix1d zz = result;
    var thisElems = _elements;
    var yElems = yy._elements;
    var resultElems = zz._elements;
    if (thisElems == null || yElems == null || resultElems == null) {
      throw 'Null elements should be impossible.';
    }
    int thiss = _columnStride;
    int ys = yy.stride;
    int results = zz.stride;

    int indexThis = index(0, 0);
    int indexY = yy.index(0);
    int indexResult = zz.index(0);

    int cols = columns;
    for (int row = rows; row >= 0; row--) {
      double sum = 0;

      // loop unrolled
      int i = indexThis - thiss;
      int j = indexY - ys;
      for (int k = (cols % 4) - 1; k >= 0; k--) {
        sum += thisElems[i += thiss] * yElems[j += ys];
      }
      for (int k = (cols / 4).floor() - 1; k >= 0; k--) {
        sum += thisElems[i += thiss] * yElems[j += ys] +
          thisElems[i += thiss] * yElems[j += ys] +
          thisElems[i += thiss] * yElems[j += ys] +
          thisElems[i += thiss] * yElems[j += ys];
      }

      resultElems[indexResult] = alpha * sum + beta * resultElems[indexResult];
      indexThis += _rowStride;
      indexResult += results;
    }

    return result;
  }

  /**
   * Linear algebraic matrix-matrix multiplication; result = alpha * this x b + beta * result.
   */
  Matrix2d mul2dWithParams(Matrix2d b, Matrix2d result, double alpha, double beta,
      bool transposeThis, bool transposeB) {
    if (transposeThis) {
      return viewTranspose().mul2dWithParams(b, result, alpha, beta, false, transposeB);
    }
    /*if (b is SparseMatrix2d || b is RCMatrix2d) {
      if (result == null) {
        return b.mul2dWithParams(this, null, alpha, beta, !transposeB, true).viewTranspose();
      } else {
        b.mul2dWithParams(this, result.viewTranspose(), alpha, beta, !transposeB, true);
        return result;
      }
    }*/

    if (transposeB) {
      return mul2dWithParams(b.viewTranspose(), result, alpha, beta, transposeThis, false);
    }

    int m = rows;
    int n = columns;
    int p = b.columns;
    if (result == null) {
      result = new DenseMatrix2d.withDimensions(m, p);
    }
    if (result is! DenseMatrix2d) {
      return super.mul2dWithParams(b, result, alpha, beta,transposeThis, transposeB);
    }
    if (b.rows != n) {
      throw new ArgumentError('Matrix2D inner dimensions must agree:${toStringShort()}, ' +
          b.toStringShort());
    }
    if (result.rows != m || result.columns != p) {
      throw new ArgumentError('Incompatible result matrix: ${toStringShort()}, ' +
          '${b.toStringShort()}, ${result.toStringShort()}');
    }
    if (this == result || b == result) {
      throw new ArgumentError('Matrices must not be identical');
    }

    DenseMatrix2d bb = b as DenseMatrix2d;
    DenseMatrix2d cc = result as DenseMatrix2d;
    var thisElems = this._elements;
    var bElems = bb._elements;
    var resultElems = result._elements;
    if (thisElems == null || bElems == null || resultElems == null) {
      throw 'Null elements should be impossible.';
    }

    int cThis = _columnStride;
    int cB = bb._columnStride;
    int cResult = cc._columnStride;

    int rThis = _rowStride;
    int rB = bb._rowStride;
    int rResult = cc._rowStride;

    /*
    A is blocked to hide memory latency
        xxxxxxx B
        xxxxxxx
        xxxxxxx
    A
    xxx     xxxxxxx C
    xxx     xxxxxxx
    ---     -------
    xxx     xxxxxxx
    xxx     xxxxxxx
    ---     -------
    xxx     xxxxxxx
    */
    int m_optimal = ((blockSize - n) / (n + 1)).floor();
    if (m_optimal <= 0) {
      m_optimal = 1;
    }
    int blocks = (m / m_optimal).floor();
    int rr = 0;
    if (m % m_optimal != 0) {
      blocks++;
    }
    while (--blocks >= 0) {
      int jB = bb.index(0, 0);
      int indexA = index(rr,0);
      int jC =  cc.index(rr,0);
      rr += m_optimal;
      if (blocks == 0) {
        m_optimal += m - rr;
      }

      for (int j = p - 1; j >= 0; j--) {
        int iA = indexA;
        int iC = jC;
        for (int i = m_optimal - 1; i >= 0; i--) {
          int kA = iA;
          int kB = jB;
          double s = 0.0;

          /*
          // Not unrolled:
          for (int k = n; --k >= 0; ) {
            //s += getQuick(i,k) * B.getQuick(k,j);
            s += AElems[kA] * BElems[kB];
            kB += rB;
            kA += cA;
          }
          */

          // Loop unrolled:
          kA -= cThis;
          kB -= rB;

          for (int k = (n % 4) - 1; k >= 0; k--) {
            s += thisElems[kA += cThis] * bElems[kB += rB];
          }
          for (int k = (n / 4).floor() - 1; k >= 0; k--) {
            s += thisElems[kA += cThis] * bElems[kB += rB] +
              thisElems[kA += cThis] * bElems[kB += rB] +
              thisElems[kA += cThis] * bElems[kB += rB] +
              thisElems[kA += cThis] * bElems[kB += rB];
          }

          resultElems[iC] = alpha * s + beta * resultElems[iC];
          iA += rThis;
          iC += rResult;
        }
        jB += cB;
        jC += cResult;
      }
    }
    return result;
  }

  /**
   * Returns the sum of all cells.
   */
  double sum() {
    double sum = 0;
    if (_elements == null) {
      throw 'Null elements should be impossible.';
    }
    int idx = index(0, 0);
    for (int row = rows - 1; row >= 0; row--) {
      for (int i = idx, column = columns; column >= 0; column--) {
        sum += _elements[i];
        i += _columnStride;
      }
      idx += _rowStride;
    }
    return sum;
  }
}
