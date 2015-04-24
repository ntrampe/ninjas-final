//
//  Filename:     matrix_triangular_base.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                matrix_triangular_base class
//


template <class T>
matrix_triangular_base<T>::matrix_triangular_base()
{
	this->m_size = 1;

	this->initMatrix();
}


template <class T>
matrix_triangular_base<T>::matrix_triangular_base(const size_t aSize)
{
	this->m_size = aSize;

	this->initMatrix();
}


template <class T>
matrix_triangular_base<T>::matrix_triangular_base(const matrix_triangular_base<T>& aCopy)
{
	*this = aCopy;
}


template <class T>
matrix_triangular_base<T>::~matrix_triangular_base()
{
	this->m_size = 0;
}


template <class T>
size_t matrix_triangular_base<T>::rows() const
{
	return m_size;
}


template <class T>
size_t matrix_triangular_base<T>::columns() const
{
	return m_size;
}


template <class T>
size_t matrix_triangular_base<T>::memorySize() const
{
	return ((1 + m_size) * m_size) / 2;
}


template <class T>
size_t matrix_triangular_base<T>::size() const
{
	return m_size;
}


template <class T>
matrix_triangular_base<T>& matrix_triangular_base<T>::operator=(const matrix_triangular_base<T>& aRHS)
{
	matrix_base<T>::operator=(aRHS);
	return *this;
}


template <class U>
vector<U> operator*(const vector<U>& aLHS, const matrix_triangular_base<U>& aRHS)
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
vector<U> operator*(const matrix_triangular_base<U>& aLHS, const vector<U>& aRHS)
{
	return (aRHS * aLHS);
}


template <class T>
T& matrix_triangular_base<T>::at(const size_t aRow, const size_t aColumn)
{
	if (!withinTriangle(aRow, aColumn))
	{
		return m_outsideElement;
	}

	return matrix_base<T>::at(aRow, aColumn);
}


template <class T>
T matrix_triangular_base<T>::at(const size_t aRow, const size_t aColumn) const
{
	if (!withinTriangle(aRow, aColumn))
	{
		return m_outsideElement;
	}

	return matrix_base<T>::at(aRow, aColumn);
}


template <class T>
void matrix_triangular_base<T>::setupMatrix(const size_t aRow, const size_t aColumn)
{
	m_size = std::min(aRow, aColumn);
}


template <class T>
void matrix_triangular_base<T>::initMatrix(const size_t aMemorySize)
{
	matrix_base<T>::initMatrix(aMemorySize);

	m_outsideElement = T();
}

