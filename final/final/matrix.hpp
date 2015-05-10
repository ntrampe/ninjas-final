//
//  Filename:     matrix.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the matrix class
//
//


template <class T>
matrix<T>::matrix()
{
	m_rows = m_columns = 1;

	this->initMatrix();
}


template <class T>
matrix<T>::matrix(const size_t aRows, const size_t aColumns)
{
	m_rows = aRows;
	m_columns = aColumns;

	this->initMatrix();
}


template <class T>
matrix<T>::matrix(const matrix<T>& aCopy)
{
	*this = aCopy;
}


template <class T>
matrix<T>::matrix(const matrix_base<T>& aMatrix)
{
	m_rows = aMatrix.rows();
	m_columns = aMatrix.columns();

	this->copyMatrix(aMatrix);
}


template <class T>
matrix<T>::matrix(const vector<T>& aVector)
{
	m_rows = aVector.size();
	m_columns = 1;

	this->initMatrix();

	for (size_t i = 0; i < m_rows; i++)
	{
		this->at(i, 0) = aVector[(unsigned int)i];
	}
}


template <class T>
matrix<T>::~matrix()
{
	m_rows = 0;
	m_columns = 0;
}


template <class T>
size_t matrix<T>::rows() const
{
	return m_rows;
}


template <class T>
size_t matrix<T>::columns() const
{
	return m_columns;
}


template <class T>
size_t matrix<T>::memorySize() const
{
	return m_rows * m_columns;
}


template <class T>
std::string matrix<T>::name() const
{
	return "Dense Matrix";
}


template <class T>
void matrix<T>::insertVectorAtRow(const vector<T>& aVector, const size_t aRow)
{
	if (aVector.size() != m_columns)
	{
		throw std::invalid_argument("matrix: provided vector size == matrix columns");
	}

	// resize matrix to make room for aVector
	if (aRow > m_rows)
	{
		this->resize(m_rows + (aRow - m_rows) + 1, m_columns);
	}
	else
	{
		this->resize(m_rows + 1, m_columns);
	}


	// shift rows down as necessary
	for (size_t i = m_rows-1; i > aRow; i--)
	{
		for (size_t j = 0; j < m_columns; j++)
		{
			this->at(i, j) = this->at(i-1, j);
		}
	}

	// set elements of inserted vector
	for (size_t i = 0; i < m_columns; i++)
	{
		this->at(aRow, i) = aVector[i];
	}
}


template <class T>
void matrix<T>::insertVectorAtColumn(const vector<T>& aVector, const size_t aColumn)
{
	if (aVector.size() != m_rows)
	{
		throw std::invalid_argument("matrix: provided vector size == matrix rows");
	}

	// resize matrix to make new room for aVector
	if (aColumn > m_columns)
	{
		this->resize(m_rows, m_columns + (aColumn - m_columns) + 1);
	}
	else
	{
		this->resize(m_rows, m_columns + 1);
	}


	// shift columns over as necessary
	for (size_t i = m_columns-1; i > aColumn; i--)
	{
		for (size_t j = 0; j < m_rows; j++)
		{
			this->at(j, i) = this->at(j, i-1);
		}
	}

	// set elements of inserted vector
	for (size_t i = 0; i < m_rows; i++)
	{
		this->at(i, aColumn) = aVector[i];
	}
}

template <class T>
matrix<T> matrix<T>::getSubMatrix(const size_t aStartRow, const size_t aStartColumn, const size_t aEndRow, const size_t aEndColumn) const
{
	if (aStartRow > m_rows || aStartColumn > m_columns || aEndRow > m_rows || aEndColumn > m_columns)
	{
		throw std::invalid_argument("matrix: dimensions of sub matrix must be within dimensions of parent matrix");
	}

	if (aStartRow == aEndRow || aStartColumn == aEndColumn)
	{
		throw std::invalid_argument("matrix: dimensions of sub matrix must not be zero");
	}

	size_t startRow = fmin(aStartRow, aEndRow);
	size_t endRow = fmax(aStartRow, aEndRow);
	size_t startColumn = fmin(aStartColumn, aEndColumn);
	size_t endColumn = fmax(aStartColumn, aEndColumn);

	matrix<T> res(endRow - startRow, endColumn - startColumn);

	for (size_t i = startRow; i < endRow; i++)
	{
		for (size_t j = startColumn; j < endColumn; j++)
		{
			res(i - startRow, j - startColumn) = this->at(i, j);
		}
	}

	return res;
}


template <class T>
double matrix<T>::determinant() const
{
	if (m_rows != m_columns)
	{
		throw std::invalid_argument("matrix: Matrix must be a square matrix in order to calculate the determinant");
	}

	double res = 0;
	bool pos = true;
	size_t col = 0;
	matrix<T>* sub;

	if (m_rows == 1)
	{
		res = this->at(0,0);
	}
	else
	{
		sub = new matrix<T>(m_rows-1,m_rows-1);

		for (size_t i = 0; i < m_rows; i++)
		{
			col = 0;

			for (size_t j = 0; j < m_rows; j++)
			{
				if (j != i)
				{
					for (size_t k = 1; k < m_rows; k++)
					{
						(*sub)(k-1, col) = this->at(k, j);
					}
					col++;
				}
			}

			res += (pos == true ? 1 : -1) * this->at(0,i) * sub->determinant();

			pos = !pos;
		}

		delete sub;
	}

	return res;
}


template <class T>
matrix<T>& matrix<T>::operator=(const matrix_base<T>& aRHS)
{
	matrix_base<T>::operator=(aRHS);
	return *this;
}


template <class U>
matrix<U> operator-(const matrix_base<U>& aRHS)
{
	matrix<U> res = aRHS;

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.m_data[i] = -res.m_data[i];
	}

	return res;
}


template <class U>
matrix<U> operator~(const matrix_base<U>& aRHS)
{
	matrix<U> res(aRHS.columns(), aRHS.rows());

	size_t col = 0;
	size_t row = 0;

	for (size_t i = 0; i < res.rows() * res.columns(); i++)
	{
		res.convertIndexToCoordinates(row, col, i);

		res(row, col) = aRHS(col, row);
	}

	return res;
}


template <class U>
matrix<U> operator+(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS)
{
	if (aLHS.rows() != aRHS.rows() || aLHS.columns() != aRHS.columns())
	{
		throw std::invalid_argument("matrix: Matrix addition must be performed with same size matrices");
	}

	matrix<U> res = aLHS;
	size_t col = 0;
	size_t row = 0;

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.convertIndexToCoordinates(row, col, i);

		res(row, col) += aRHS(row, col);
	}

	return res;
}


template <class U>
matrix<U> operator-(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS)
{
	return (aLHS + -aRHS);
}


template <class U>
matrix<U> operator*(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS)
{
	if (aLHS.rows() != aRHS.columns() || aLHS.columns() != aRHS.rows())
	{
		throw std::invalid_argument("matrix: Matrix multiplication must be performed with nAm x mBn matricies");
	}

	matrix<U> res(aLHS.rows(), aLHS.rows());
	size_t col = 0;
	size_t row = 0;

	for (size_t i = 0; i < res.memorySize(); i++)
	{
		res.convertIndexToCoordinates(row, col, i);

		for (size_t j = 0; j < aLHS.columns(); j++)
		{
			res(row, col) += aLHS(row, j) * aRHS(j, col);
		}
	}

	return res;
}


template <class U>
matrix<U> operator*(const double& aLHS, const matrix_base<U>& aRHS)
{
	matrix<U> res = aRHS;

	size_t col = 0;
	size_t row = 0;

	for (size_t i = 0; i < res.rows() * res.columns(); i++)
	{
		res.convertIndexToCoordinates(row, col, i);

		res(row, col) *= aLHS;
	}

	return res;
}


template <class U>
matrix<U> operator*(const matrix_base<U>& aLHS, const double& aRHS)
{
	return aRHS * aLHS;
}


template <class U>
vector<U> operator*(const vector<U>& aLHS, const matrix_base<U>& aRHS)
{
	if (aLHS.size() != aRHS.columns())
	{
		throw std::invalid_argument("matrix: Matrix-vector multiplication must be performed where vector size equals matrix columns");
	}

	vector<U> res;

	for (size_t i = 0; i < aRHS.rows(); i++)
	{
		res.push_back(U());
	}

	for (size_t i = 0; i < aRHS.rows(); i++)
	{
		for (size_t j = 0; j < aRHS.columns(); j++)
		{
			res[i] += aLHS[j] * aRHS(i, j);
		}
	}

	return res;
}


template <class U>
vector<U> operator*(const matrix_base<U>& aLHS, const vector<U>& aRHS)
{
	return aRHS * aLHS;
}


template <class T>
void matrix<T>::convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const
{
	if (aRow >= m_rows)
	{
		throw std::out_of_range("matrix: provided rows >= matrix rows");
	}

	if (aColumn >= m_columns)
	{
		throw std::out_of_range("matrix: provided columns >= matrix columns");
	}

	aIndex = (aRow * m_columns + aColumn);
}


template <class T>
void matrix<T>::convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const
{
	if (aIndex >= memorySize())
	{
		throw std::out_of_range("matrix: provided index >= matrix data array size");
	}

	aColumn = aIndex % m_columns;
	aRow = aIndex / m_columns;
}


template <class T>
void matrix<T>::setupMatrix(const size_t aRow, const size_t aColumn)
{
	m_rows = aRow;
	m_columns = aColumn;
}
