//
//  Filename:     matrix_tridiagonal.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the matrix_tridiagonal class
//
//


template<class T>
matrix_tridiagonal<T>::matrix_tridiagonal(const matrix_tridiagonal<T>& aCopy) : matrix_banded<T>(aCopy)
{
	*this = aCopy;
}


template<class T>
matrix_tridiagonal<T>::matrix_tridiagonal(const matrix_base<T>& aMatrix)
{
	this->copyMatrix(aMatrix);
}


template<class T>
std::string matrix_tridiagonal<T>::name() const
{
	return "Tridiagonal Matrix";
}


template <class T>
bool matrix_tridiagonal<T>::solveMatrix(const vector<T>& aB, vector<T>& aX)
{
	//// Thomas Algorithm

	vector<T> a, b, c, d;

	b = aB;
	aX.clear();
	aX.reserve(this->columns(), true);

	for (size_t i = 0; i < this->columns(); i++)
	{
		d.push_back(this->at(i,i));
	}

	for (size_t i = 0; i < this->rows()-1; i++)
	{
		c.push_back(this->at(i, i+1));
	}

	for (size_t i = 1; i < this->rows(); i++)
	{
		a.push_back(this->at(i, i-1));
	}

	for (size_t i = 1; i < this->columns(); i++)
	{
		d[i] = d[i] - (a[i - 1] / d[i - 1])*c[i - 1];
		b[i] = b[i] - (a[i - 1] / d[i - 1])*b[i - 1];
	}

	aX[this->columns()-1] = b[this->columns() - 1] / d[this->columns() - 1];

	for (int i = static_cast<int>(this->columns()) - 2; i >= 0; i--)
	{
		aX[i] = (b[i] - c[i]*aX[i+1]) / d[i];
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


template<class T>
matrix_tridiagonal<T>& matrix_tridiagonal<T>::operator=(const matrix_base<T>& aRHS)
{
	matrix_base<T>::operator=(aRHS);
	return *this;
}

