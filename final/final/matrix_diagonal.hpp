//
//  Filename:     matrix_diagonal.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the matrix_diagonal class
//
//


template<class T>
matrix_diagonal<T>::matrix_diagonal(const matrix_diagonal<T>& aCopy) : matrix_banded<T>(aCopy)
{
	*this = aCopy;
}


template<class T>
matrix_diagonal<T>::matrix_diagonal(const matrix_base<T>& aMatrix)
{
	this->copyMatrix(aMatrix);
}


template<class T>
std::string matrix_diagonal<T>::name() const
{
	return "Diagonal Matrix";
}


template <class T>
void matrix_diagonal<T>::replaceVectorAtRow(const vector<T>& aVector, const size_t aRow)
{
	if (aRow >= this->m_size)
	{
		throw std::out_of_range("matrix: provided row >= matrix rows");
	}

	if (aVector.size() != this->m_size)
	{
		throw std::invalid_argument("matrix: provided vector size == matrix columns");
	}

	this->m_data[aRow] = aVector[aRow];
}


template <class T>
void matrix_diagonal<T>::replaceVectorAtColumn(const vector<T>& aVector, const size_t aColumn)
{
	if (aColumn >= this->m_size)
	{
		throw std::out_of_range("matrix: provided column >= matrix column");
	}

	if (aVector.size() != this->m_size)
	{
		throw std::invalid_argument("matrix: provided vector size == matrix rows");
	}

	this->m_data[aColumn] = aVector[aColumn];
}


template<class T>
matrix_diagonal<T>& matrix_diagonal<T>::operator=(const matrix_base<T>& aRHS)
{
	matrix_base<T>::operator=(aRHS);
	return *this;
}


template<class U>
matrix_diagonal<U> operator*(const matrix_diagonal<U>& aLHS, const matrix_diagonal<U>& aRHS)
{
	if (aLHS.rows() != aRHS.columns() || aLHS.columns() != aRHS.rows())
	{
		throw std::invalid_argument("matrix_diagonal: Matrix multiplication must be performed with nAm x mBn matricies");
	}

	matrix_diagonal<U> res = aLHS;

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] *= aRHS.m_data[i];
	}

	return res;
}

template<class U>
vector<U> operator*(const vector<U>& aLHS, const matrix_diagonal<U>& aRHS)
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
vector<U> operator*(const matrix_diagonal<U>& aLHS, const vector<U>& aRHS)
{
	return aRHS * aLHS;
}



template<class T>
void matrix_diagonal<T>::convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const
{
	if (aRow != aColumn)
	{
		throw std::out_of_range("matrix_diagonal: provided row != provided column");
	}

	aIndex = aRow;
}


template<class T>
void matrix_diagonal<T>::convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const
{
	if (aIndex >= this->memorySize())
	{
		throw std::out_of_range("matrix_diagonal: provided index >= memory size");
	}

	aRow = aIndex;
	aColumn = aIndex;
}

