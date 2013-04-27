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


/// An aggregator function for matrix aggregation.
typedef double Aggregator(double total, double immediate);


/**
 * Abstract base class for 1-d matrices.
 * A matrix has a number of cells, which are assigned upon instance construction.
 * Elements are accessed via zero based indexes.  Legal indexes are of the form [0..size - 1].
 */
abstract class Matrix1d {

  /**
   * Returns the [double] value in the matrix cell value at coordinate [index].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  double operator[](int index);

  /**
   * Sets the [double] value in the matrix cell value at coordinate [index].  Behavior when
   * accessing elements out of bounds of the matrix is undefined.
   */
  void operator[]=(int index, double value);

  /**
   * Construct and returns a new 2-d matrix of the same dynamic type as this matrix, with specified
   * [rows] and [columns] sizes.
   * TODO(klochek): This is waiting on the Matrix2d implementation.
   */
  /*Matrix2d*/ like2D(int rows, int columns);

  /**
   * Construct and returns a new selection view.
   */
  Matrix1d viewSelectionLike(List<int> offsets);

  /// The number of cells this matrix (or view) has.
  int _size;

  /// The index of the first element.
  int _zero;

  /// The number of indexes between any two elements, i.e. index(i+1) - index(i).
  int _stride;

  /// Flag indicating whether or not this matrix is a view of another matrix.
  bool _isView = false;

  /**
   * Returns the [int] position of the given absolute rank within the (virtual or non-virtual)
   * internal 1-dimensional array.
   */
  int _offset(int absRank) {
    return absRank;
  }

  /**
   * Returns the [int] absolute rank of the given relative [rank].
   */
  int _rank(int rank) {
    return _zero + rank * _stride;
  }

  /**
   * Throws a [RangeError] if [index] is outside the bounds of the matrix.
   */
  void _checkValidIndex(int index) {
    if (index < 0 || index >= size) {
      throw new RangeError("Attempted to access ${toStringShort()} at index=$index");
    }
  }

  /**
   * Checks whether [indexes] are legal and throws an [RangeError], if necessary.
   */
  void _checkValidIndexes(List<int> indexes) {
    indexes.forEach(_checkValidIndex);
  }

  /**
   * Checks whether the matrix contains the given range and throws a [RangeError], if necessary.
   */
  void _checkValidRange(int index, int width) {
    if (index < 0 || index + width > size) {
      throw new RangeError("index: $index, width: $width, size=$size");
    }
  }

  /**
   * Checks if [B] (either a [List<double>] or a [Matrix1d] is of the same size as this matrix.
   * Throws a [RangeError] if the sizes do not match.
   */
  void _checkValidSize(dynamic B) {
    if (B is List<double>) {
      if (size != B.length) {
        throw new ArgumentError("Incompatible sizes: ${toStringShort()} and ${B.length}");
      }
    } else if (B is Matrix1d) {
      if (size != B.size) {
        throw new ArgumentError("Incompatible sizes: ${toStringShort()} and ${B.toStringShort()}");
      }
    }
  }

  /**
   * Sets up a matrix with [size] cells, and an optional [zero] index and [stride].  [size] must be
   * greater than zero.
   */
  void _setUp(int size, [int zero = 0, int stride = 1]) {
    if (size < 0) {
      throw new ArgumentError("negative size");
    }

    _size = size;
    _zero = zero;
    _stride = stride;
    _isView = false;
  }

  /**
   * Returns the number of cells.
   */
  int get size {
    return _size;
  }

  /**
   * Returns the position of the element with the given relative [rank] within this matrix.
   */
  int index(int rank) {
    return _offset(_rank(rank));
  }

  /**
   * Returns the stride of the given [dimension] (axis, rank).  [dimension] must be > 0.
   */
  int stride(int dimension) {
    if (dimension != 0) {
      throw new ArgumentError("invalid dimension: $dimension used to access ${toStringShort()}");
    }
    return _stride;
  }

  /**
   * Returns the number of cells having non-zero values, up to an optional maximum of
   * [maxCardinality]; ignores tolerance.
   */
  int cardinality([int maxCardinality = -1]) {
    int cardinality = 0;
    maxCardinality = maxCardinality > -1 ? maxCardinality : _size;
    int i = _size - 1;
    while (i >= 0 && cardinality < maxCardinality) {
      if (this[i] != 0) {
        cardinality++;
      }
    }
    return cardinality;
  }

  /**
   * Returns a string representation of the receiver's shape.
   */
  String toStringShort() {
    throw "Not implemented yet!";
    // TODO(klochek): return AbstractFormatter.shape(this);
  }

  /**
   * Executes [procedure] on every cell in the matrix.  Returns [true] if [procedure] returns [true]
   * for every single cell.
   */
  bool every(bool procedure(double x)) {
    for (int i = size; i >= 0; i--) {
      if (!procedure(this[i])) {
        return false;
      }
    }
    return true;
  }

  /**
   * Applies a function [f] to each cell and aggregates the results, via the supplied [aggr]
   * function.  Returns NAN if [size] = 0.
   */
  double aggregate(Aggregator aggr, double f(double x)) {
    if (size == 0) {
      return double.NAN;
    }
    double a = f(this[_size - 1]);
    for (int i = _size - 2; i >= 0; i--) {
      a = aggr(a, f(this[i]));
    }
    return a;
  }

  /**
   * Applies a function [f] to each corresponding cell of this matrix and [other], and aggregates
   * the results via the supplied [aggr] function.  Returns NAN if [size] = 0;
   */
  double aggregateOther(Matrix1d other, Aggregator aggr, double f(double x, double y)) {
    _checkValidSize(other);
    if (_size == 0) {
      return double.NAN;
    }
    double a = f(this[_size - 1], other[_size - 1]);
    for (int i = _size - 2; i >= 0; i--) {
      a = aggr(a, f(this[i], other[i]));
    }
    return a;
  }

  /**
   * Returns the sum of all cells.
   */
  double sum() {
    if (_size == 0) {
      return 0.0;
    }
    return aggregate((a, b) => a + b, (a) => a);
  }

  /**
   * Returns the dot product this matrix and [y];
   */
  double dot(Matrix1d y) {
    return dotRange(y, 0, size);
  }

  /**
   * Returns the dot product of this matrix and [y], only on the range of cells starting at [from]
   * and continuing for [length] cells.
   */
  double dotRange(Matrix1d y, int from, int length) {
    if (from < 0 || length <= 0) {
      return 0.0;
    }

    int tail = from + length;
    if (_size < tail) {
      tail = _size;
    }
    if (y._size < tail) {
      tail = y._size;
    }
    length = tail - from;

    double sum = 0.0;
    int i = tail - 1;
    for (int k = length - 1; k >= 0; k--, i--) {
      sum += this[i] * y[i];
    }
    return sum;
  }

  /**
   * Returns the dot product of this matrix and [y], only on the range of cells starting at [from]
   * and continuing for [length] cells.  Only non-zero values in [nonZeroIndexes] are considered.
   */
  double dotRangeWithIndices(Matrix1d y, int from, int length, List<int> nonZeroIndexes) {
    // determine minimum length
    if (from < 0 || length <= 0) {
      return 0.0;
    }

    int tail = from + length;
    if (_size < tail) {
      tail = _size;
    }
    if (y._size < tail) {
      tail = y._size;
    }
    length = tail - from - 1;
    if (length <= 0) {
      return 0.0;
    }

    // setup
    List<int> nonZeroIndexElements = nonZeroIndexes;
    int index = 0;
    int s = nonZeroIndexes.length;

    // skip to start
    while ((index < s) && nonZeroIndexElements[index] < from) {
      index++;
    }

    // now the sparse dot product
    int i;
    double sum = 0.0;
    while ((length-- >= 0) && (index < s) && ((i=nonZeroIndexElements[index]) < tail)) {
      sum += this[i] * y[i];
      index++;
    }

    return sum;
  }

  /**
   * Returns the dot product of this matrix and [y], where only the non-zero values in
   * [nonZeroIndexes] are considered.
   */
  double dotWithIndices(Matrix1d y, List<int> nonZeroIndexes) {
    return dotRangeWithIndices(y, 0, _size, nonZeroIndexes);
  }

  /**
   * Sets all cells to the state specified by [v].
   * [v] can be a [List<double>], a [Matrix1d], a [double], or a [Function].
   * If [v] is a [List<double>], then [this] gets every value in [v].  The sizes of [v] and [this]
   * must be equal.
   * If [v] is a [Matrix1d], then [this] gets every value in [v].  The sizes of [v] and [this]
   * must be equal.
   * If [v] is a [Function], then each value in [this] is transformed by [v].
   * If [v] is a [double], then every value in [this] is set to [v].
   */
  void assign(dynamic v) {
    if (v is List<double>) {
      if (v.length != _size) {
        throw new ArgumentError("Must have same number of cells: length=${v.length} size()=$_size");
      }
      for (int i = _size - 1; i >= 0; i--) {
        this[i] = v[i];
      }
    } else if (v is double) {
      for (int i = _size - 1; i >= 0; i--) {
        this[i] = v;
      }
    } else if (v is Function) {
      for (int i = _size - 1; i >= 0; i--) {
        this[i] = v(this[i]);
      }
    } else if (v is Matrix1d) {
      if (v == this) {
        return;
      }
      _checkValidSize(v);
      if (haveSharedCells(v)) {
        v = v.copy();
      }
      for (int i = _size - 1; i >= 0; i--) {
        this[i] = v[i];
      }
    }
  }

  /**
   * Assigns the result of a function [f] to each cell in [this], so that this[i] =
   * f(this[i], y[i]).
   */
  void assignWithFunc(Matrix1d y, double f(double thisVal, double thatVal)) {
    _checkValidSize(y);
    for (int i = _size - 1; i >= 0; i--) {
      this[i] = f(this[i], y[i]);
    }
  }

  /**
   * Assigns the result of a function [f] to each cell in [this], so that this[i] =
   * f(this[i], y[i]).  Uses the indices in [nonZeroIndexes] to potentially optimize the application
   * process.
   */
  void assignWithIndices(Matrix1d y, function, List<int> nonZeroIndexes) {
    _checkValidSize(y);

    // TODO(klochek): colt uses non-zero indices as a hint for certain functions (+, *).  Provide
    // some const functions that do the same, and then specialize this function for them.
    assign(y, function);
  }

  /**
   * Construct and returns a new empty matrix of the same dynamic type as this matrix, with an
   * optionally-different specified [size].
   */
  Matrix1d like([int size = - 1]) {
    size = size > -1 ? size : _size;
    return like(size);
  }

  /**
   * Swaps each element this[i] with other[i].
   */
  void swap(Matrix1d other) {
    _checkValidSize(other);
    for (int i = size - 1; i >= 0; i--) {
      double tmp = this[i];
      this[i] = other[i];
      other[i] = tmp;
    }
  }

  /**
   * Constructs and returns a deep copy of this matrix.
   */
  Matrix1d copy() {
    Matrix1d copy = like();
    copy.assign(this);
    return copy;
  }

  /**
   * Returns whether all cells are equal to the given [value].
   * [value] may be a [double], or a [Matrix1d].
   */
  bool equals(dynamic value) {
     if (value is double) {
       //return cern.colt.matrix.linalg.Property.DEFAULT.equals(this,value);
       throw "Not implemented yet!";
     }
     if (this == value) {
       return true;
     }
     if (value == null) {
       return false;
     }
     if (!(value is Matrix1d)) {
       return false;
     }
     //return cern.colt.matrix.linalg.Property.DEFAULT.equals(this,obj);
     throw "Not implemented yet!";
  }

  /**
   * Safely returns the matrix cell value at coordinate [index].  Throws a [RangeError] if the
   * supplied [index] is out of bounds.
   */
  double safeGet(int index) {
    _checkValidIndex(index);
    return this[index];
  }

  /**
   * Safely sets the matrix cell at coordinate [index] to the specified [value].  Throws a
   * [RangeError] if no such index exists.
   */
  void safeSet(int index, double value) {
    _checkValidIndex(index);
    this[index] = value;
  }

  /**
   * Returns the content of this matrix if it is a wrapper; or <tt>this</tt> otherwise.
   * Override this method in wrappers.
   */
  Matrix1d get content {
    return this;
  }

  /**
   * Fills the coordinates and values of cells having non-zero values into the specified [indexList]
   * and [valueList], up to an optional [maxCardinality].
   * In general, fill order is unspecified.
   */
  void getNonZeros(List<int> indexList, List<double> valueList, [int maxCardinality = -1]) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    maxCardinality = maxCardinality > -1 ? maxCardinality : _size;

    int card = cardinality(maxCardinality);
    if (fillIndexList) {
      indexList.length = card;
      indexList.fillRange(0, card, 0);
    }
    if (fillValueList) {
      valueList.length = card;
      valueList.fillRange(0, card, 0);
    }
    if (card >= maxCardinality) {
      return;
    }

    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    for (int i=0; i < _size; i++) {
      double value = this[i];
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }
  }

  /**
   * Returns [true] if both matrices share at least one identical cell.
   */
  bool haveSharedCells(Matrix1d other) {
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
  bool haveSharedCellsRaw(Matrix1d other) {
    // TODO(klochek): huh?  What's going on here?  Throw an exception?
    return false;
  }

  /**
   * Constructs and returns a new view equal to the receiver.
   * The returned matrix is backed by this matrix, so changes in the returned matrix are reflected
   * in this matrix, and vice-versa.
   * Use [copy()] to construct an independent deep copy rather than a new view.
   */
   Matrix1d view() {
     throw "Not implemented yet!";
   }

  /**
   * Constructs and returns a new flipped view.
   */
  Matrix1d viewFlip() {
    return view().vFlip() as Matrix1d;
  }

  /**
   * Constructs and returns a new view of this matrix that is a [width] sub matrix starting at
   * [index].
   */
  Matrix1d viewPart(int index, int width) {
    return view().vPart(index, width) as Matrix1d;
  }

  /**
   * Constructs and returns a new selection view that is a matrix holding the supplied [indexes].
   */
  Matrix1d viewSelectionOf(List<int> indexes) {
    if (indexes == null) {
      indexes = new List<int>(_size);
      for (int i = _size - 1; i >= 0; i--) {
        indexes[i] = i;
      }
    }

    List<int> offsets = new List<int>(indexes.length);
    for (int i = indexes.length - 1; i >= 0; i--) {
      offsets[i] = index(indexes[i]);
    }
    return viewSelectionLike(offsets);
  }

  /**
   * Constructs and returns a new selection view that is a matrix holding the cells matching the given
   * [condition] function.
   */
  Matrix1d viewSelectionWith(bool condition(double x)) {
    List<int> matches = new List<int>();
    for (int i = 0; i < _size; i++) {
      if (condition(this[i])) {
        matches.add(i);
      }
    }
    // TODO(klochek): huh? was: matches.trimToSize();
    return viewSelectionOf(matches);
  }

  /**
   * Sorts the vector into ascending order, according to the natural ordering.
   * This sort is guaranteed to be stable.
   */
  Matrix1d viewSorted() {
    throw "Not implemented yet!";
    //return cern.colt.matrix.doublealgo.Sorting.mergeSort.sort(this);
  }

  /**
   * Constructs and returns a new stride view which is a sub matrix consisting of every [stride]-th cell.
   */
  Matrix1d viewStrides(int stride) {
    return view().vStrides(stride) as Matrix1d;
  }

  /**
   * Self-modifying version of [viewFlip()]. What used to be index 0 is now index size - 1, and
   * vice-versa.
   */
  void vFlip() {
    if (_size > 0) {
      _zero += (_size - 1) * _stride;
      _stride = -_stride;
      _isView = true;
    }
  }

  /**
   * Self-modifying version of [viewPart()].
   */
  void vPart(int index, int width) {
    _checkValidRange(index, width);
    _zero += _stride * index;
    _size = width;
    _isView = true;
  }

  /**
   * Self modifying version of [viewStrides()].
   * @throws IndexOutOfBoundsException if <tt>stride <= 0</tt>.
   */
  void vStrides(int stride) {
    if (stride <= 0) {
      throw new RangeError("illegal stride: $stride");
    }
    _stride *= stride;
    if (_size != 0) {
      _size = ((_size - 1) / stride + 1).floor();
    }
    _isView = true;
  }

  /**
   * Fills the cell values in the optionally specified list [values], returning it.  If a list is
   * specified, it must be of the same size as this matrix.
   */
  List<double> toArray([List<double> values]) {
    if (values != null) {
      if (values.length < size) {
        throw new ArgumentError("values too small");
      }
    } else {
      values = new List<double>(size);
    }
    for (int i = size - 1; i >= 0; i--) {
      values[i] = this[i];
    }
    return values;
  }

  /**
   * Returns a string representation using default formatting.
   */
  String toString() {
    throw "Not implemented yet!";
    //return new cern.colt.matrix.doublealgo.Formatter().toString(this);
  }
}