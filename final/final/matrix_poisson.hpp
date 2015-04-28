//
//  Filename:     matrix_poisson.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                matrix_poisson class
//

#include "matrix_triangular_lower.h"

template <class T>
matrix_poisson<T>::matrix_poisson()
{
	setupMatrix(1, 1);

	this->initMatrix();
}


template <class T>
matrix_poisson<T>::matrix_poisson(const size_t aSlices)
{
  size_t s = 0;
  m_slices = aSlices;
  s = size();
  
	setupMatrix(s, s);

	this->initMatrix();
}


template <class T>
matrix_poisson<T>::matrix_poisson(const matrix_poisson<T>& aCopy)
{
	*this = aCopy;
}


template <class T>
matrix_poisson<T>::~matrix_poisson()
{
	m_slices = 0;
}


template <class T>
size_t matrix_poisson<T>::rows() const
{
	return size();
}


template <class T>
size_t matrix_poisson<T>::columns() const
{
	return size();
}


template <class T>
size_t matrix_poisson<T>::memorySize() const
{
  return  matrix_base<T>::lengthOfDiagonal(0) + 
          matrix_base<T>::lengthOfDiagonal(1) - (m_slices-1-1) +
          matrix_base<T>::lengthOfDiagonal(m_slices-1);
}


template <class T>
size_t matrix_poisson<T>::size() const
{
	return (m_slices - 1)*(m_slices - 1);
}


template <class T>
size_t matrix_poisson<T>::slices() const
{
  return m_slices;
}


template <class T>
size_t matrix_poisson<T>::meshSize() const
{
  return (m_slices - 1)*(m_slices - 1);
}


template <class T>
size_t matrix_poisson<T>::band() const
{
  return (m_slices - 1);
}


template <class T>
std::string matrix_poisson<T>::name() const
{
	return "Poisson Matrix";
}


template <class T>
void matrix_poisson<T>::replaceVectorAtRow(const vector<T>& aVector, const size_t aRow)
{
	if (aRow >= this->size())
	{
		throw std::out_of_range("matrix: provided row >= matrix rows");
	}

	if (aVector.size() != this->size())
	{
		throw std::invalid_argument("matrix: provided vector size == matrix columns");
	}

	for (size_t i = 0; i <= aRow; i++)
	{
		this->at(aRow, i) = aVector[i];
	}
}


template <class T>
void matrix_poisson<T>::replaceVectorAtColumn(const vector<T>& aVector, const size_t aColumn)
{
	if (aColumn >= this->size())
	{
		throw std::out_of_range("matrix: provided column >= matrix column");
	}

	if (aVector.size() != this->size())
	{
		throw std::invalid_argument("matrix: provided vector size == matrix rows");
	}

	for (size_t i = 0; i <= aColumn; i++)
	{
		this->at(aColumn, i) = aVector[i];
	}
}


template <class T>
bool matrix_poisson<T>::solveMatrix(const vector<T>& aB, vector<T>& aX)
{
	//// Cholesky Decomposition

	matrix_poisson<T> aL(*this);
	matrix_poisson<T> aA;
	vector<T> y;
	double sum = 0;

	y.reserve(size(), true);
	y = aB;

	for (int k = 0; k < static_cast<int>(size()); k++)
	{
		for (int i = 0; i <= k-1; i++)
		{
			sum = 0;
			for (int j = 0; j <= i-1; j++)
			{
				sum = sum + aL(i,j) * aL(k,j);
			}
			aL(k,i) = (aL(k,i) - sum)/aL(i,i);
		}

		sum = 0;

		for (int j = 0; j <= k-1; j++)
		{
			sum = sum + aL(k,j) * aL(k,j);
		}

		aL(k,k) = sqrt(aL(k,k) - sum);
	}

	aA = aL;

	for (size_t col = 0; col < aA.size(); col++)
	{
		if (aA(col, col) != 0 && aA(col, col) != 1)
		{
			y[col] = y[col] / aA(col, col);

			aA(col, col) = 1;
		}

		for (size_t row = col + 1; row < aA.size(); row++)
		{
			y[row] = -aA(row, col) * y[col] + y[row];

			aA(row, col) = 0;
		}
	}

	double csum = 0;

	// create solution array
	aX.clear();
	aX.reserve(aL.columns(), true);

	// backward substitution
	// solve for the unknown variables from the bottom up
	for (int row = (int)aL.columns() - 1; row >= 0; row--)
	{
		csum = 0;

		// compute sum of the coefficients times the old
		// x values
		for (size_t col = row; col < aL.columns(); col++)
		{
			csum += aL(row, col) * aX[col];
		}

		// compute the new x by subtracting the sum from the
		// right-hand side of the equation and dividing by the
		// coefficient corresponding to the current x
		//
		// xi = (an - (a1 + a2 + ... + an-1 (not including ai))) / ai
		//
		// where i = row
		aX[row] = (y[row] - csum) / aL(row, row);

		// fix rounding errors (e.g. 0 != 0.000000000001)
		// for display purposes
		if (equivalent(aX[row], 0))
		{
			aX[row] = 0;
		}
	}

	// check for no solution
	// could move this check into the back substitution loop,
	// but this is prettier
	for (size_t i = 0; i < aX.size(); i++)
	{
		if (aX[i] != aX[i])
		{
			// aX[i] is not a number
			// no solution
			return false;
		}
	}

	// aX now contains the solution
	return true;
}


template <class T>
matrix_poisson<T>& matrix_poisson<T>::operator=(const matrix_base<T>& aRHS)
{
	matrix_base<T>::operator=(aRHS);
	return *this;
}


template <class U>
matrix_poisson<U> operator-(const matrix_poisson<U>& aRHS)
{
	matrix_poisson<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = -aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_poisson<U> operator+(const matrix_poisson<U>& aLHS, const matrix_poisson<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix_poisson: matrices must be the same size to perform matrix addition");
	}

	matrix_poisson<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = aLHS.m_data[i] + aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_poisson<U> operator-(const matrix_poisson<U>& aLHS, const matrix_poisson<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix_poisson: matrices must be the same size to perform matrix subtraction");
	}

	matrix_poisson<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = aLHS.m_data[i] - aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_poisson<U> operator*(const double& aLHS, const matrix_poisson<U>& aRHS)
{
	matrix_poisson<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = aLHS * aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_poisson<U> operator*(const matrix_poisson<U>& aLHS, const double& aRHS)
{
	return (aRHS * aLHS);
}


template <class U>
vector<U> operator*(const vector<U>& aLHS, const matrix_poisson<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix_poisson: Matrix-vector multiplication must be performed where vector size equals matrix columns");
	}

	vector<U> res;

	for (size_t i = 0; i < aRHS.size(); i++)
	{
		res.push_back(U());
	}

	for (size_t i = 0; i < aRHS.rows(); i++)
	{
		for (size_t j = 0; j <= i; j++)
		{
			res[i] += aLHS[j] * aRHS(i, j);
		}
	}

	return res;
}


template <class U>
vector<U> operator*(const matrix_poisson<U>& aLHS, const vector<U>& aRHS)
{
	return (aRHS * aLHS);
}



template <class T>
void matrix_poisson<T>::convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const
{
	if (!this->withinDiagonal(aRow, aColumn))
	{
		throw std::out_of_range("matrix_poisson: provided row, column must be within diagonal");
	}

	size_t realRow = aRow;
	size_t realCol = aColumn;

	if (aRow < aColumn)
	{
		realRow = aColumn;
		realCol = aRow;
	}

	aIndex = 0;

//  for (size_t i = 0; i < realRow; i++)
//  {
//    aIndex += widthAtRow(i);
//  }
//  
//  for (size_t i = startAtRow(realRow); i < realCol; i++)
//  {
//    aIndex++;
//  }
  
  
  if (realRow == realCol)
  {
    aIndex = realRow;
  }
  else if (realRow == realCol+1 && (realRow % (m_slices - 1) != 0))
  {
    aIndex = matrix_base<T>::lengthOfDiagonal(0);
    
    for (size_t i = 1; i < realRow; i++)
    {
      if (i % (m_slices - 1) != 0)
      {
        aIndex++;
      }
    }
  }
  else if (realRow == realCol+(m_slices - 1))
  {
    aIndex =  matrix_base<T>::lengthOfDiagonal(0) + 
              matrix_base<T>::lengthOfDiagonal(1) - (m_slices-1-1);
    
    for (size_t i = (m_slices - 1); i < realRow; i++)
    {
      aIndex++;
    }
  }
}


template <class T>
void matrix_poisson<T>::convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const
{
	if (aIndex >= this->memorySize())
	{
		throw std::out_of_range("matrix_poisson: provided index >= memory size");
	}

	aRow = 0;
	aColumn = 0;

  for (size_t i = 0; i < aIndex; i++)
  {
    aColumn++;
    
    if (aColumn == endAtRow(aRow)+1)
    {
      aRow++;
      aColumn = startAtRow(aRow);
    }
  }
}


template<class T>
bool matrix_poisson<T>::withinDiagonal(const size_t aRow, const size_t aColumn) const
{
  if (!this->withinDimensions(aRow, aColumn))
  {
    return false;
  }
  
  size_t realRow = aRow;
  size_t realCol = aColumn;
  
  if (aRow < aColumn)
  {
    realRow = aColumn;
    realCol = aRow;
  }
  
  if ((realRow != realCol) && !((realRow == realCol+1 && (realRow % (m_slices - 1) != 0))) && !(realRow == realCol+(m_slices - 1)))
  {
    return false;
  }
  
  return true;
}


template <class T>
T& matrix_poisson<T>::at(const size_t aRow, const size_t aColumn)
{
  if (!withinDiagonal(aRow, aColumn))
  {
    return m_outsideElement;
  }
  
	if (aRow < aColumn)
	{
		return matrix_base<T>::at(aColumn, aRow);
	}

	return matrix_base<T>::at(aRow, aColumn);
}


template <class T>
T matrix_poisson<T>::at(const size_t aRow, const size_t aColumn) const
{
  if (!withinDiagonal(aRow, aColumn))
  {
    return m_outsideElement;
  }
  
	if (aRow < aColumn)
	{
		return matrix_base<T>::at(aColumn, aRow);
	}

	return matrix_base<T>::at(aRow, aColumn);
}


template <class T>
void matrix_poisson<T>::setupMatrix(const size_t aRow, const size_t aColumn)
{
  m_outsideElement = T();
}


template <class T>
void matrix_poisson<T>::initMatrix(const size_t aMemorySize)
{
  matrix_base<T>::initMatrix(aMemorySize);
  size_t s = size();
  
  for (size_t i = 0; i < memorySize(); i++)
  {
    if (i < s)
    {
      this->m_data[i] = 4;
    }
    else
    {
      this->m_data[i] = -1;
    }
  }
  
  m_outsideElement = T();
}


template<class T>
size_t matrix_poisson<T>::bandWidth() const
{
  return band() + 1;
}


template<class T>
size_t matrix_poisson<T>::widthAtRow(const size_t aRow) const
{
  size_t res = 0;
  
  if (aRow < band())
  {
    res = bandWidth() - band() + aRow;
  }
  else
  {
    res = bandWidth();
  }
  
  return res;
}


template<class T>
size_t matrix_poisson<T>::startAtRow(const size_t aRow) const
{
  size_t res = 0;
  
  if (aRow <= band())
  {
    res = 0;
  }
  else
  {
    res = aRow - band();
  }
  
  return res;
}


template<class T>
size_t matrix_poisson<T>::endAtRow(const size_t aRow) const
{
  size_t res = 0;
  
  if (aRow > (size() - 1) - band())
  {
    res = size() - 1;
  }
  else
  {
    res = aRow + band();
  }
  
  return res;
}


template<class T>
size_t matrix_poisson<T>::widthAtColumn(const size_t aColumn) const
{
  return widthAtRow(aColumn);
}


template<class T>
size_t matrix_poisson<T>::startAtColumn(const size_t aColumn) const
{
  return startAtRow(aColumn);
}


template<class T>
size_t matrix_poisson<T>::endAtColumn(const size_t aColumn) const
{
  return endAtRow(aColumn);
}


