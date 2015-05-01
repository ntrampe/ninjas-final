//
//  Filename:     vector.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the vector class
//
//


template <class T>
vector<T>::vector(const vector<T>& aCopy)
{
	*this = aCopy;
}


template <class T>
vector<T>::~vector()
{
	if (m_data != NULL)
	{
		delete [] m_data;
		m_data = NULL;
	}
}


template <class T>
void vector<T>::clear()
{
	if (m_data != NULL)
	{
		m_size = 0;
		m_maxSize = 0;
		delete [] m_data;
		m_data = NULL;
	}
}


template <class T>
void vector<T>::reserve(const size_t aQuantity, const bool aUpdateSize)
{
	setMemory(aQuantity);

	if (aUpdateSize)
	{
		m_size = 0;

		for (size_t i = 0; i < aQuantity; i++)
		{
			push_back(T());
		}
	}
}


template <class T>
void vector<T>::push_back(const T& aElement)
{
	addMemory();
	m_data[m_size] = aElement;
	m_size++;
}


template <class T>
void vector<T>::remove(const T& aElement)
{
	int i = 0;
	bool found = false;

	while (i < m_size && !found)
	{
		if (m_data[i] == aElement)
			found = true;
		else i++;
	}

	if (found)
	{
		for (int j = i; j < m_size - 1; j++)
			m_data[j] = m_data[j + 1];
		m_size--;
	}
}


template <class T>
void vector<T>::pop_back()
{
	if (m_size > 0)
	{
		m_size--;
		removeMemory();
	}
	else
	{
		throw std::length_error("empty vector");
	}
}


template <class T>
size_t vector<T>::size() const
{
	return m_size;
}


template <class T>
size_t vector<T>::maxSize() const
{
	return m_maxSize;
}


template <class T>
bool vector<T>::isMultiple(const vector<T>& aRHS) const
{
	if (m_size != aRHS.m_size)
		return false;

	//TODO: find a better way to accomplish this...

	T max, min;
	int mult = std::numeric_limits<int>::max();

	for (unsigned int i = 0; i < m_size; i++)
	{
		max = std::max(aRHS[i], m_data[i]);
		min = std::min(aRHS[i], m_data[i]);

		if (fmod(max, min) != 0)
		{
			// the elements are not multiples of each other
			return false;
		}
		else
		{
			if (mult == std::numeric_limits<int>::max())
			{
				// first multiple catch, set the multiple
				if (min != 0)
				{
					mult = max / min;
				}
			}
			else
			{
				if (min != 0)
				{
					if (mult != max / min)
					{
						// elements were multiples, but not the same multiple
						return false;
					}
				}
			}
		}
	}

	// all elements are the same multiple of each other
	return true;
}


template <class T>
T& vector<T>::operator[](size_t aIndex)
{
	if (aIndex >= m_size)
		throw std::out_of_range("index >= m_size");
	return m_data[aIndex];
}


template <class T>
const T& vector<T>::operator[](size_t aIndex) const
{
	if (aIndex >= m_size)
		throw std::out_of_range("index >= m_size");
	return m_data[aIndex];
}


template <class U>
std::ostream& operator<<(std::ostream& aOutput, const vector<U>& aVector)
{
	std::stringstream ss;

	ss << "< ";

	for (unsigned int i = 0; i < aVector.size(); i++)
	{
		ss << aVector[i] << (i != aVector.size() - 1 ? ", " : "");
	}

	ss << " >";

	aOutput << ss.str();

	return aOutput;
}


template <class U>
std::istream& operator>>(std::istream& aInput, vector<U>& aVector)
{
	for (size_t i = 0; i < aVector.size(); i++)
	{
		aInput >> aVector[i];
	}

	return aInput;
}


template <class U>
vector<U> operator+(const vector<U>& aLHS, const vector<U>& aRHS)
{
	vector<U> res;
	unsigned int size = static_cast<unsigned int>(std::max(aLHS.size(), aRHS.size()));

	res.reserve(size);

	for (unsigned int i = 0; i < size; i++)
	{
		//ternary checks make sure we don't exceed the bounds of either array
		//if we exceed the bounds of one array, just use a default U object
		res.push_back((i < aLHS.size() ? aLHS[i] : U()) + (i < aRHS.size() ? aRHS[i] : U()));
	}

	return res;
}


template <class U>
vector<U> operator-(const vector<U>& aLHS, const vector<U>& aRHS)
{
	return aLHS + -aRHS;
}


template <class U>
U operator*(const vector<U>& aLHS, const vector<U>& aRHS)
{
	U res = U();
	unsigned int size = static_cast<unsigned int>(std::max(aLHS.size(), aRHS.size()));

	for (unsigned int i = 0; i < size; i++)
	{
		//ternary checks make sure we don't exceed the bounds of either array
		//if we exceed the bounds of one array, just use a default U object
		res += ((i < aLHS.size() ? aLHS[i] : U()) * (i < aRHS.size() ? aRHS[i] : U()));
	}

	return res;
}


template <class U>
vector<U> operator*(const double aLHS, const vector<U>& aRHS)
{
	vector<U> res;

	for (unsigned int i = 0; i < aRHS.size(); i++)
	{
		res.push_back(aLHS * aRHS[i]);
	}

	return res;
}

template <class U>
vector<U> operator*(const vector<U>& aLHS, const double aRHS)
{
	return aRHS * aLHS;
}


template <class U>
vector<U> operator-(const vector<U>& aVector)
{
	vector<U> res;

	res.reserve(aVector.size());

	for (unsigned int i = 0; i < aVector.size(); i++)
	{
		res.push_back(-aVector[i]);
	}

	return res;
}


template <class U>
bool operator==(const vector<U>& aLHS, const vector<U>& aRHS)
{
	if (aLHS.size() != aRHS.size())
		return false;

	for (unsigned int i = 0; i < aLHS.size(); i++)
	{
		if (aLHS[i] != aRHS[i])
		{
			return false;
		}
	}

	return true;
}


template <class U>
bool operator!=(const vector<U>& aLHS, const vector<U>& aRHS)
{
	return !(aLHS == aRHS);
}


template <class T>
vector<T>& vector<T>::operator=(const vector<T>& aRHS)
{
	if (this != &aRHS)
	{
		m_size = aRHS.m_size;
		m_maxSize = aRHS.m_maxSize;
		m_data = new T [m_maxSize];
		for(unsigned int k = 0; k < m_maxSize; k++)
			m_data[k] = aRHS.m_data[k];
	}
	return *this;
}


template <class T>
void vector<T>::addMemory()
{
	if (m_size+1 == m_maxSize)
	{
		setMemory(m_maxSize*2);
	}
	else if (m_size == 0 && m_data == NULL)
	{
		m_maxSize = 2;
		m_data = new T [m_maxSize];
	}
}


template <class T>
void vector<T>::removeMemory()
{
	if (static_cast<float>(m_size)/m_maxSize <= 0.25 && m_size != 0)
	{
		setMemory(m_maxSize/2);
	}
	else if (m_size == 0)
	{
		clear();
	}
}


template <class T>
void vector<T>::setMemory(const size_t aMemory)
{
	unsigned int mem = 1;

	while (mem < aMemory)
		mem *= 2;

	if (m_data != NULL)
	{
		T * tempData = new T [m_maxSize];

		for (unsigned int x = 0; x < m_size; x++)
			tempData[x] = m_data[x];

		m_maxSize = mem;

		delete [] m_data;
		m_data = new T [m_maxSize];

		for (unsigned int x = 0; x < m_size; x++)
			m_data[x] = tempData[x];

		delete [] tempData;
	}
	else
	{
		m_maxSize = mem;
		m_data = new T [m_maxSize];
	}
}
