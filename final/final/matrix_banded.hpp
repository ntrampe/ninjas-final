//
//  Filename:     matrix_banded.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the matrix_banded class
//
//

#include "matrix.h"

template<class T>
matrix_banded<T>::matrix_banded()
{
	m_band = 0;
	setupMatrix(1, 1);

	this->initMatrix();
}


template<class T>
matrix_banded<T>::matrix_banded(const size_t aSize, const size_t aBand)
{
	m_band = aBand;

	setupMatrix(aSize, aSize);

	this->initMatrix();
}


template<class T>
matrix_banded<T>::matrix_banded(const matrix_banded<T>& aCopy)
{
	*this = aCopy;
}


template<class T>
matrix_banded<T>::matrix_banded(const matrix_base<T>& aMatrix)
{
	this->copyMatrix(aMatrix);
}


template <class T>
matrix_banded<T>::~matrix_banded()
{
	m_size = 0;
}


template<class T>
size_t matrix_banded<T>::rows() const
{
	return m_size;
}


template<class T>
size_t matrix_banded<T>::columns() const
{
	return m_size;
}


template<class T>
size_t matrix_banded<T>::size() const
{
	return m_size;
}


template<class T>
size_t matrix_banded<T>::band() const
{
	return m_band;
}


template<class T>
size_t matrix_banded<T>::memorySize() const
{
	size_t res = 0;

	for (int i = -static_cast<int>(m_band); i <= static_cast<int>(m_band); i++)
	{
		res += m_size - abs(i);
	}

	return res;
}


template<class T>
std::string matrix_banded<T>::name() const
{
	return "Banded Matrix";
}


template<class T>
std::string matrix_banded<T>::description() const
{
	std::stringstream ss;

	ss << ", band(" << band() << ")";

	return matrix_base<T>::description().append(ss.str());
}


template <class T>
void matrix_banded<T>::replaceVectorAtRow(const vector<T>& aVector, const size_t aRow)
{
	if (aRow >= m_size)
	{
		throw std::out_of_range("matrix: provided row >= matrix rows");
	}

	if (aVector.size() != m_size)
	{
		throw std::invalid_argument("matrix: provided vector size == matrix columns");
	}

	for (size_t i = startAtRow(aRow); i <= endAtRow(aRow); i++)
	{
		this->at(aRow, i) = aVector[i];
	}
}


template <class T>
void matrix_banded<T>::replaceVectorAtColumn(const vector<T>& aVector, const size_t aColumn)
{
	if (aColumn >= m_size)
	{
		throw std::out_of_range("matrix: provided column >= matrix column");
	}

	if (aVector.size() != m_size)
	{
		throw std::invalid_argument("matrix: provided vector size == matrix rows");
	}

	for (size_t i = startAtColumn(aColumn); i <= endAtColumn(aColumn); i++)
	{
		this->at(i, aColumn) = aVector[i];
	}
}


template <class T>
bool matrix_banded<T>::solveMatrix(const vector<T>& aB, vector<T>& aX)
{
	//ehhhh

	matrix<T> m(*this);

	return m.solveMatrix(aB, aX);
}


template<class T>
matrix_banded<T>& matrix_banded<T>::operator=(const matrix_banded<T>& aRHS)
{
	m_band = aRHS.band();
	matrix_base<T>::operator=(aRHS);
	return *this;
}


template<class U>
matrix_banded<U> operator-(const matrix_banded<U>& aRHS)
{
	matrix_banded<U> res = aRHS;

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = -res.m_data[i];
	}

	return res;
}


template<class U>
matrix_banded<U> operator+(const matrix_banded<U>& aLHS, const matrix_banded<U>& aRHS)
{
	if (aLHS.rows() != aRHS.rows() && aLHS.columns() != aRHS.columns())
	{
		throw std::invalid_argument("matrix_banded: matrices must be the same size to perform matrix addition");
	}

	matrix_banded<U> res = aLHS;

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] += aRHS.m_data[i];
	}

	return res;
}


template<class U>
matrix_banded<U> operator-(const matrix_banded<U>& aLHS, const matrix_banded<U>& aRHS)
{
	if (aLHS.rows() != aRHS.rows() && aLHS.columns() != aRHS.columns())
	{
		throw std::invalid_argument("matrix_banded: matrices must be the same size to perform matrix subtraction");
	}

	matrix_banded<U> res = aLHS;

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] -= aRHS.m_data[i];
	}

	return res;
}


template<class U>
matrix_banded<U> operator*(const double& aLHS, const matrix_banded<U>& aRHS)
{
	matrix_banded<U> res = aRHS;

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] *= aLHS;
	}

	return res;
}


template<class U>
matrix_banded<U> operator*(const matrix_banded<U>& aLHS, const double& aRHS)
{
	return aRHS * aLHS;
}


template<class U>
vector<U> operator*(const vector<U>& aLHS, const matrix_banded<U>& aRHS)
{
	if (aLHS.size() != aRHS.columns())
	{
		throw std::invalid_argument("matrix: Matrix-vector multiplication must be performed where vector size equals matrix columns");
	}

	vector<U> res;

	for (size_t i = 0; i < aRHS.size(); i++)
	{
		res.push_back(aLHS[i] * aRHS.m_data[i]);
	}

	return res;
}


template<class U>
vector<U> operator*(const matrix_banded<U>& aLHS, const vector<U>& aRHS)
{
	return aRHS * aLHS;
}



template<class T>
void matrix_banded<T>::convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const
{
	if (!withinDiagonal(aRow, aColumn))
	{
		throw std::out_of_range("matrix_banded: provided row, column must be within the diagonal");
	}

	aIndex = 0;

	for (size_t i = 0; i < aRow; i++)
	{
		aIndex += widthAtRow(i);
	}

	for (size_t i = startAtRow(aRow); i < aColumn; i++)
	{
		aIndex++;
	}
}


template<class T>
void matrix_banded<T>::convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const
{
	if (aIndex >= this->memorySize())
	{
		throw std::out_of_range("matrix_banded: provided index >= memory size");
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
bool matrix_banded<T>::withinDiagonal(const size_t aRow, const size_t aColumn) const
{
	if (!this->withinDimensions(aRow, aColumn))
	{
		return false;
	}

	return (aRow <= aColumn + m_band && aColumn <= aRow + m_band);
}


template<class T>
T& matrix_banded<T>::at(const size_t aRow, const size_t aColumn)
{
	if (!withinDiagonal(aRow, aColumn))
	{
		return m_outsideElement;
	}

	return matrix_base<T>::at(aRow, aColumn);
}


template<class T>
T matrix_banded<T>::at(const size_t aRow, const size_t aColumn) const
{
	if (!withinDiagonal(aRow, aColumn))
	{
		return m_outsideElement;
	}

	return matrix_base<T>::at(aRow, aColumn);
}


template<class T>
void matrix_banded<T>::setupMatrix(const size_t aRow, const size_t aColumn)
{
	if (aRow < m_band + 1 || aColumn < m_band + 1)
	{
		throw std::out_of_range("matrix_banded: size must be >= bandwidth");
	}

	m_size = std::max(aRow, aColumn);
	m_outsideElement = T();
}


template<class T>
size_t matrix_banded<T>::bandWidth() const
{
	return 2 * m_band + 1;
}


template<class T>
size_t matrix_banded<T>::widthAtRow(const size_t aRow) const
{
	size_t res = 0;

	if (aRow < m_band)
	{
		res = bandWidth() - m_band + aRow;
	}
	else if (aRow > (m_size - 1) - m_band)
	{
		res = bandWidth() - (m_band - ((m_size - 1) - aRow));
	}
	else
	{
		res = bandWidth();
	}

	return res;
}


template<class T>
size_t matrix_banded<T>::startAtRow(const size_t aRow) const
{
	size_t res = 0;

	if (aRow <= m_band)
	{
		res = 0;
	}
	else
	{
		res = aRow - m_band;
	}

	return res;
}


template<class T>
size_t matrix_banded<T>::endAtRow(const size_t aRow) const
{
	size_t res = 0;

	if (aRow > (m_size - 1) - m_band)
	{
		res = m_size - 1;
	}
	else
	{
		res = aRow + m_band;
	}

	return res;
}


template<class T>
size_t matrix_banded<T>::widthAtColumn(const size_t aColumn) const
{
	return widthAtRow(aColumn);
}


template<class T>
size_t matrix_banded<T>::startAtColumn(const size_t aColumn) const
{
	return startAtRow(aColumn);
}


template<class T>
size_t matrix_banded<T>::endAtColumn(const size_t aColumn) const
{
	return endAtRow(aColumn);
}



