//
//  Filename:     matrix_triangular_upper.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                matrix_triangular_upper class
//


template <class T>
std::string matrix_triangular_upper<T>::name() const
{
	return "Upper Triangular Matrix";
}


template <class T>
void matrix_triangular_upper<T>::replaceVectorAtRow(const vector<T>& aVector, const size_t aRow)
{
	if (aRow >= this->m_size)
	{
		throw std::out_of_range("matrix: provided row >= matrix rows");
	}

	if (aVector.size() != this->m_size)
	{
		throw std::invalid_argument("matrix: provided vector size == matrix columns");
	}

	for (size_t i = aRow; i < this->m_size; i++)
	{
		this->at(aRow, i) = aVector[i];
	}
}


template <class T>
void matrix_triangular_upper<T>::replaceVectorAtColumn(const vector<T>& aVector, const size_t aColumn)
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
		this->at(i, aColumn) = aVector[i];
	}
}


template <class T>
matrix_triangular_upper<T>& matrix_triangular_upper<T>::operator=(const matrix_base<T>& aRHS)
{
	matrix_base<T>::operator=(aRHS);
	return *this;
}


template <class U>
matrix_triangular_upper<U> operator-(const matrix_triangular_upper<U>& aRHS)
{
	matrix_triangular_upper<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = -aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_triangular_upper<U> operator+(const matrix_triangular_upper<U>& aLHS, const matrix_triangular_upper<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix_triangular_upper: matrices must be the same size to perform matrix addition");
	}

	matrix_triangular_upper<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = aLHS.m_data[i] + aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_triangular_upper<U> operator-(const matrix_triangular_upper<U>& aLHS, const matrix_triangular_upper<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix_triangular_upper: matrices must be the same size to perform matrix addition");
	}

	matrix_triangular_upper<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = aLHS.m_data[i] - aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_triangular_upper<U> operator*(const double& aLHS, const matrix_triangular_upper<U>& aRHS)
{
	matrix_triangular_upper<U> res(aRHS.size());

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = aLHS * aRHS.m_data[i];
	}

	return res;
}


template <class U>
matrix_triangular_upper<U> operator*(const matrix_triangular_upper<U>& aLHS, const double& aRHS)
{
	return (aRHS * aLHS);
}


template <class U>
matrix_triangular_upper<U> operator*(const matrix_triangular_upper<U>& aLHS, const matrix_triangular_upper<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix: Matrix multiplication must be performed with nAm x mBn matricies");
	}

	matrix_triangular_upper<U> res(aLHS.size());

	for (size_t row = 0; row < res.size(); row++)
	{
		for (size_t col = row; col < res.size(); col++)
		{
			for (size_t i = row; i <= col; i++)
			{
				res(row, col) += aLHS(row, i) * aRHS(i, col);
			}
		}
	}

	return res;
}


template <class U>
vector<U> operator*(const vector<U>& aLHS, const matrix_triangular_upper<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
	{
		throw std::invalid_argument("matrix: Matrix-vector multiplication must be performed where vector size equals matrix columns");
	}

	vector<U> res;

	for (size_t i = 0; i < aRHS.size(); i++)
	{
		res.push_back(U());
	}

	for (size_t i = 0; i < aRHS.rows(); i++)
	{
		for (size_t j = i; j < aRHS.columns(); j++)
		{
			res[i] += aLHS[j] * aRHS(i, j);
		}
	}

	return res;
}


template <class U>
vector<U> operator*(const matrix_triangular_upper<U>& aLHS, const vector<U>& aRHS)
{
	return (aRHS * aLHS);
}


template <class T>
void matrix_triangular_upper<T>::convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const
{
	if (!withinTriangle(aRow, aColumn))
	{
		throw std::out_of_range("matrix: provided row <= provided column");
	}

	aIndex = 0;

	for (size_t i = 0; i < aRow; i++)
	{
		aIndex += (this->m_size - i - 1);
	}

	for (size_t i = 0; i < aColumn; i++)
	{
		aIndex++;
	}
}


template <class T>
void matrix_triangular_upper<T>::convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const
{
	if (aIndex >= this->memorySize())
	{
		throw std::out_of_range("matrix: provided index >= memory size");
	}

	aRow = 0;
	aColumn = 0;

	for (size_t i = 0; i < aIndex; i++)
	{
		aColumn++;

		if (aColumn == this->m_size)
		{
			aRow++;
			aColumn = aRow;
		}
	}
}


template <class T>
bool matrix_triangular_upper<T>::withinTriangle(const size_t aRow, const size_t aColumn) const
{
	if (!this->withinDimensions(aRow, aColumn))
		return false;

	return (aRow <= aColumn);
}
