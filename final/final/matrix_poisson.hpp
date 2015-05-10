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
  return 2;
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
void matrix_poisson<T>::resize(const size_t aRows, const size_t aColumns)
{
  if (aRows != aColumns)
  {
    throw std::invalid_argument("matrix_poisson: provided rows != columns");
  }
  
  if (sqrt(aRows) != floor(sqrt(aRows)))
  {
    
  }
  
  m_slices = sqrt(aRows) + 1;
}


template <class T>
void matrix_poisson<T>::replaceVectorAtRow(const vector<T>& aVector, const size_t aRow)
{
	if (aRow >= this->size())
	{
		throw std::out_of_range("matrix_poisson: provided row >= matrix rows");
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
matrix_poisson<T>& matrix_poisson<T>::operator=(const matrix_poisson<T>& aRHS)
{
  m_slices = aRHS.slices();
	matrix_base<T>::operator=(aRHS);
	return *this;
}


template <class T>
void matrix_poisson<T>::convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const
{
	if (!this->withinData(aRow, aColumn))
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
  
  if (realRow == realCol)
  {
    aIndex = 0;
  }
  else
  {
    aIndex = 1;
  }
}


template <class T>
void matrix_poisson<T>::convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const
{
  aRow = 0;
  aColumn = 0;
  
  throw std::logic_error("matrix_poisson: call to convertIndexToCoordinates is ambiguous as it's not a one-to-one function");
}


template<class T>
bool matrix_poisson<T>::withinData(const size_t aRow, const size_t aColumn) const
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
  if (!withinData(aRow, aColumn))
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
  if (!withinData(aRow, aColumn))
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
  
  this->m_data[0] = T(1);
  this->m_data[1] = T(-0.25);
  
  m_outsideElement = T();
}


