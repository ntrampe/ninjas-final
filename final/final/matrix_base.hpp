//
//  Filename:     matrix_base.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                matrix_base class
//


template <class T>
matrix_base<T>::~matrix_base()
{
	if (m_data != nullptr)
	{
		delete [] m_data;
		m_data = nullptr;
	}
}


template <class T>
std::string matrix_base<T>::description() const
{
	std::stringstream ss;

	ss << name() << " - dim(" << rows() << " x " << columns() << "), size(" << memorySize() << ")";

	return ss.str();
}


template <class T>
void matrix_base<T>::initMatrix(const size_t aMemorySize)
{
	if (m_data != nullptr)
	{
		delete [] m_data;
		m_data = nullptr;
	}

	size_t mem = (aMemorySize == 0 ? memorySize() : aMemorySize);

	m_data = new T[mem];

	for (size_t i = 0; i < mem; i++)
	{
		m_data[i] = T();
	}
}


template <class T>
void matrix_base<T>::copyMatrix(const matrix_base<T>& aMatrix)
{
	setupMatrix(aMatrix.rows(), aMatrix.columns());
	initMatrix();

	for (size_t row = 0; row < rows(); row++)
	{
		for (size_t col = 0; col < columns(); col++)
		{
			at(row, col) = aMatrix(row, col);
		}
	}
}


template <class T>
size_t matrix_base<T>::longestElement() const
{
	size_t res = 0;
	std::stringstream ss;

	for (size_t i = 0; i < memorySize(); i++)
	{
		ss.clear();
		ss.str(std::string());
		ss << m_data[i];

		if (ss.str().size() > res)
		{
			res = ss.str().size();
		}
	}

	return res;
}


template <class T>
void matrix_base<T>::clear(const T aDefault)
{
	for (size_t i = 0; i < memorySize(); i++)
	{
		m_data[i] = aDefault;
	}
}


template <class T>
void matrix_base<T>::randomize(const unsigned int aMaximum)
{
	for (size_t i = 0; i < memorySize(); i++)
	{
		m_data[i] = static_cast<T>((rand() % aMaximum) + 1);
	}
}


template <class T>
void matrix_base<T>::order(const T aStart)
{
	for (size_t i = 0; i < memorySize(); i++)
	{
		m_data[i] = aStart + T(i);
	}
}


template <class T>
void matrix_base<T>::printMemory() const
{
	size_t mem = memorySize();
	std::cout << "[ ";

	for (size_t i = 0; i < mem; i++)
	{
		std::cout << m_data[i] << (i == mem - 1 ? " ]" : ", ");
	}

	std::cout << std::endl;
}


template <class T>
void matrix_base<T>::roundToZero()
{
	for (size_t i = 0; i < memorySize(); i++)
	{
		if (equivalent(m_data[i], T()))
		{
			m_data[i] = T();
		}
	}
}


template <class T>
void matrix_base<T>::resize(const size_t aRows, const size_t aColumns)
{
	size_t mem = memorySize();
	T* temp = new T[mem];
	point2d<size_t>* oldLoc = new point2d<size_t>[mem];
	size_t oldArea = mem;
	size_t col = 0;
	size_t row = 0;

	if (this->m_data != nullptr)
	{
		for (unsigned int i = 0; i < memorySize(); i++)
		{
			convertIndexToCoordinates(row, col, i);

			temp[i] = this->m_data[i];
			oldLoc[i].set(row, col);
		}

		delete [] this->m_data;
		this->m_data = nullptr;
	}

	this->setupMatrix(aRows, aColumns);
	this->initMatrix();

	for (size_t i = 0; i < oldArea; i++)
	{
		if (withinDimensions(oldLoc[i].x(), oldLoc[i].y()))
		{
			this->at(oldLoc[i].x(), oldLoc[i].y()) = temp[i];
		}
	}

	delete [] oldLoc;
	delete [] temp;
}


template <class T>
size_t matrix_base<T>::lengthOfDiagonal(const size_t aRow) const
{
  return columns() - aRow;
}


template <class T>
bool matrix_base<T>::isEqualTo(const matrix_base<T>& aMatrix) const
{
	if (rows() != aMatrix.rows() || columns() != aMatrix.columns())
		return false;

	for (size_t row = 0; row < aMatrix.rows(); row++)
	{
		for (size_t col = 0; col < aMatrix.columns(); col++)
		{
			if (at(row, col) != aMatrix(row, col))
			{
				return false;
			}
		}
	}

	return true;
}


template <class T>
void matrix_base<T>::vectorAtRow(const size_t aRow, vector<T>& aVector) const
{
	if (aRow >= rows())
	{
		throw std::out_of_range("matrix_base: provided row >= matrix rows");
	}

  aVector.clear();

	for (size_t i = 0; i < columns(); i++)
	{
		aVector.push_back(this->at(aRow, i));
	}
}


template <class T>
void matrix_base<T>::vectorAtColumn(const size_t aColumn, vector<T>& aVector) const
{
	if (aColumn >= columns())
	{
		throw std::out_of_range("matrix_base: provided column >= matrix column");
	}

  aVector.clear();

	for (size_t i = 0; i < rows(); i++)
	{
		aVector.push_back(this->at(i, aColumn));
	}
}


template <class T>
void matrix_base<T>::replaceVectorAtRow(const vector<T>& aVector, const size_t aRow)
{
	if (aRow >= rows())
	{
		throw std::out_of_range("matrix: provided row >= matrix rows");
	}

	if (aVector.size() != columns())
	{
		throw std::invalid_argument("matrix: provided vector size == matrix columns");
	}

	for (size_t i = 0; i < columns(); i++)
	{
		this->at(aRow, i) = aVector[i];
	}
}


template <class T>
void matrix_base<T>::replaceVectorAtColumn(const vector<T>& aVector, const size_t aColumn)
{
	if (aColumn >= columns())
	{
		throw std::out_of_range("matrix: provided column >= matrix column");
	}

	if (aVector.size() != rows())
	{
		throw std::invalid_argument("matrix: provided vector size == matrix rows");
	}

	for (size_t i = 0; i < rows(); i++)
	{
		this->at(i, aColumn) = aVector[i];
	}
}


template <class T>
bool matrix_base<T>::isRowEchelonForm() const
{
	size_t leadingCol = columns();
	size_t previousLeadingCol = columns();

	for (size_t row = 0; row < rows(); row++)
	{
		previousLeadingCol = leadingCol;
		leadingCol = columns();

		for (size_t col = 0; col < columns(); col++)
		{
			if (at(row, col) != 0)
			{
				if (leadingCol == columns())
				{
					leadingCol = col;

					if (row > 0)
					{
						if (leadingCol <= previousLeadingCol)
						{
							return false;
						}
					}
				}
			}
		}
	}

	return true;
}


template <class T>
void matrix_base<T>::switchRows(const size_t aSource, const size_t aDestination)
{
	if (aSource >= rows() || aDestination >= rows())
	{
		throw std::out_of_range("matrix_base: row switching must occur within matrix");
	}

	if (aSource != aDestination)
	{
    vector<T> temp, dest;
    this->vectorAtRow(aSource, temp);
    this->vectorAtRow(aDestination, dest);
		this->replaceVectorAtRow(dest, aSource);
		this->replaceVectorAtRow(temp, aDestination);
	}
}


template <class T>
void matrix_base<T>::scaleRow(const size_t aSource, const double aScalar, const size_t aDestination)
{
	if (aSource >= rows() || aDestination >= rows())
	{
		throw std::out_of_range("matrix_base: row scaling must occur within matrix");
	}
  
  vector<T> temp;
  this->vectorAtRow(aSource, temp);
	this->replaceVectorAtRow(aScalar * temp, aDestination);
}


template <class T>
void matrix_base<T>::addRow(const size_t aLHS, const size_t aRHS, const size_t aDestination)
{
	if (aLHS >= rows() || aRHS >= rows() || aDestination >= rows())
	{
		throw std::out_of_range("matrix_base: row addition must occur within matrix");
	}

	replaceVectorAtRow(this->vectorAtRow(aLHS) + this->vectorAtRow(aRHS), aDestination);
}


template <class T>
T& matrix_base<T>::operator()(const size_t aRow, const size_t aColumn)
{
	return at(aRow, aColumn);
}


template <class T>
T matrix_base<T>::operator()(const size_t aRow, const size_t aColumn) const
{
	return at(aRow, aColumn);
}


template <class T>
T& matrix_base<T>::operator()(const point2d<size_t>& aPoint)
{
	return at(aPoint.x(), aPoint.y());
}


template <class T>
T matrix_base<T>::operator()(const point2d<size_t>& aPoint) const
{
	return at(aPoint.x(), aPoint.y());
}


template <class T>
matrix_base<T>& matrix_base<T>::operator=(const matrix_base<T>& aRHS)
{
	if (this != &aRHS)
	{
		this->copyMatrix(aRHS);
	}

	return *this;
}


template <class U>
bool operator==(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS)
{
	if (aLHS.memorySize() != aRHS.memorySize())
		return false;

	for (size_t i = 0; i < aLHS.memorySize(); i++)
	{
		if (aLHS.m_data[i] != aRHS.m_data[i])
		{
			return false;
		}
	}

	return true;
}


template <class U>
bool operator!=(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS)
{
	return !(aLHS == aRHS);
}


template <class U>
std::ostream& operator<<(std::ostream& aOutput, const matrix_base<U>& aMatrix)
{
	int wid = std::max(4, static_cast<int>(aMatrix.longestElement()));

	for (size_t row = 0; row < aMatrix.rows(); row++)
	{
		for (size_t col = 0; col < aMatrix.columns(); col++)
		{
			aOutput << std::setw(wid) << aMatrix(row, col) << " ";
		}

		aOutput << std::endl << std::endl;
	}

	return aOutput;
}


template <class U>
std::istream& operator>>(std::istream& aInput, matrix_base<U>& aMatrix)
{
	for (unsigned int row = 0; row < aMatrix.rows(); row++)
	{
		for (unsigned int col = 0; col < aMatrix.columns(); col++)
		{
			aInput >> aMatrix(row, col);
		}
	}

	return aInput;
}


template <class T>
bool matrix_base<T>::withinDimensions(const size_t aRow, const size_t aColumn) const
{
	return (aRow < rows() && aColumn < columns());
}


template <class T>
T& matrix_base<T>::at(const size_t aRow, const size_t aColumn)
{
	if (!withinDimensions(aRow, aColumn))
	{
		throw std::out_of_range("matrix_base: provided row, column >= matrix rows,columns");
	}

	size_t idx = 0;

	convertCoordinatesToIndex(idx, aRow, aColumn);

	return m_data[idx];
}


template <class T>
T matrix_base<T>::at(const size_t aRow, const size_t aColumn) const
{
	if (!withinDimensions(aRow, aColumn))
	{
		throw std::out_of_range("matrix_base: provided row, column >= matrix rows,columns");
	}

	size_t idx = 0;

	convertCoordinatesToIndex(idx, aRow, aColumn);

	return m_data[idx];
}
