//
//  Filename:     matrix_symmetrical.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                matrix_symmetrical class
//

#include "matrix_triangular_lower.h"

template <class T>
matrix_symmetrical<T>::matrix_symmetrical()
{
	setupMatrix(1, 1);

	this->initMatrix();
}


template <class T>
matrix_symmetrical<T>::matrix_symmetrical(const size_t aSize)
{
	setupMatrix(aSize, aSize);

	this->initMatrix();
}


template <class T>
matrix_symmetrical<T>::matrix_symmetrical(const matrix_symmetrical<T>& aCopy)
{
	*this = aCopy;
}


template <class T>
matrix_symmetrical<T>::~matrix_symmetrical()
{
	m_size = 0;
}


template <class T>
size_t matrix_symmetrical<T>::rows() const
{
	return m_size;
}


template <class T>
size_t matrix_symmetrical<T>::columns() const
{
	return m_size;
}


template <class T>
size_t matrix_symmetrical<T>::memorySize() const
{
	return ((1 + m_size) * m_size) / 2;
}


template <class T>
size_t matrix_symmetrical<T>::size() const
{
	return m_size;
}


template <class T>
std::string matrix_symmetrical<T>::name() const
{
	return "Symmetrical Matrix";
}


template <class T>
void matrix_symmetrical<T>::replaceVectorAtRow(const vector<T>& aVector, const size_t aRow)
{
	if (aRow >= this->m_size)
	{
		throw std::out_of_range("matrix: provided row >= matrix rows");
	}

	if (aVector.size() != this->m_size)
	{
		throw std::invalid_argument("matrix: provided vector size == matrix columns");
	}

	for (size_t i = 0; i <= aRow; i++)
	{
		this->at(aRow, i) = aVector[i];
	}
}


template <class T>
void matrix_symmetrical<T>::replaceVectorAtColumn(const vector<T>& aVector, const size_t aColumn)
{
	if (aColumn >= this->m_size)
	{
		throw std::out_of_range("matrix: provided column >= matrix column");
	}

	if (aVector.size() != this->m_size)
	{
		throw std::invalid_argument("matrix: provided vector size == matrix rows");
	}

	for (size_t i = 0; i <= aColumn; i++)
	{
		this->at(aColumn, i) = aVector[i];
	}
}


template <class T>
bool matrix_symmetrical<T>::solveMatrix(const vector<T>& aB, vector<T>& aX)
{
	//// Cholesky Decomposition

	matrix_symmetrical<T> aL(*this);
	matrix_symmetrical<T> aA;
	vector<T> y;
	double sum = 0;

	y.reserve(m_size, true);
	y = aB;

	for (int k = 0; k < static_cast<int>(m_size); k++)
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
matrix_symmetrical<T>& matrix_symmetrical<T>::operator=(const matrix_base<T>& aRHS)
{
	matrix_base<T>::operator=(aRHS);
	return *this;
}


template <class U>
matrix_symmetrical<U> operator-(const matrix_symmetrical<U>& aRHS)
{
	matrix_symmetrical<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = -aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_symmetrical<U> operator+(const matrix_symmetrical<U>& aLHS, const matrix_symmetrical<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix_symmetrical: matrices must be the same size to perform matrix addition");
	}

	matrix_symmetrical<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = aLHS.m_data[i] + aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_symmetrical<U> operator-(const matrix_symmetrical<U>& aLHS, const matrix_symmetrical<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix_symmetrical: matrices must be the same size to perform matrix subtraction");
	}

	matrix_symmetrical<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = aLHS.m_data[i] - aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_symmetrical<U> operator*(const double& aLHS, const matrix_symmetrical<U>& aRHS)
{
	matrix_symmetrical<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = aLHS * aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_symmetrical<U> operator*(const matrix_symmetrical<U>& aLHS, const double& aRHS)
{
	return (aRHS * aLHS);
}


template <class U>
vector<U> operator*(const vector<U>& aLHS, const matrix_symmetrical<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix_symmetrical: Matrix-vector multiplication must be performed where vector size equals matrix columns");
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
vector<U> operator*(const matrix_symmetrical<U>& aLHS, const vector<U>& aRHS)
{
	return (aRHS * aLHS);
}



template <class T>
void matrix_symmetrical<T>::convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const
{
	if (!this->withinDimensions(aRow, aColumn))
	{
		throw std::out_of_range("matrix_symmetrical: provided row, column must be within matrix");
	}

	size_t realRow = aRow;
	size_t realCol = aColumn;

	if (aRow < aColumn)
	{
		realRow = aColumn;
		realCol = aRow;
	}

	aIndex = 0;

	for (size_t i = 0; i < realRow; i++)
	{
		aIndex += i+1;
	}

	for (size_t i = 0; i < realCol; i++)
	{
		aIndex++;
	}
}


template <class T>
void matrix_symmetrical<T>::convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const
{
	if (aIndex >= this->memorySize())
	{
		throw std::out_of_range("matrix_symmetrical: provided index >= memory size");
	}

	aRow = 0;
	aColumn = 0;

	for (size_t i = 0; i < aIndex; i++)
	{
		aColumn++;

		if (aColumn == (aRow + 1))
		{
			aRow++;
			aColumn = 0;
		}
	}
}


template <class T>
T& matrix_symmetrical<T>::at(const size_t aRow, const size_t aColumn)
{
	if (aRow < aColumn)
	{
		return matrix_base<T>::at(aColumn, aRow);
	}

	return matrix_base<T>::at(aRow, aColumn);
}


template <class T>
T matrix_symmetrical<T>::at(const size_t aRow, const size_t aColumn) const
{
	if (aRow < aColumn)
	{
		return matrix_base<T>::at(aColumn, aRow);
	}

	return matrix_base<T>::at(aRow, aColumn);
}


template <class T>
void matrix_symmetrical<T>::setupMatrix(const size_t aRow, const size_t aColumn)
{
	m_size = std::max(aRow, aColumn);
}


