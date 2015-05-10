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


template<class T>
matrix_tridiagonal<T>& matrix_tridiagonal<T>::operator=(const matrix_base<T>& aRHS)
{
	matrix_base<T>::operator=(aRHS);
	return *this;
}

