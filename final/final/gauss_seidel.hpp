//
//  Filename:     gauss_seidel.hpp
//  Programmer:   James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   7 - A Parameterized Gaussian-Seidel Iteration Class
//
//  Description:  Here are the definitions of the functions declared in the gauss_seidel class
//
//

template <class T>
gauss_seidel<T>::gauss_seidel()
{
	m_error_tol = 0.0000001;
}

template <class T>
gauss_seidel<T>::~gauss_seidel()
{

}

template <class T>
gauss_seidel<T>::gauss_seidel(const double aErrorTolerance)
{
	m_error_tol = aErrorTolerance;
}


template <class T>
bool gauss_seidel<T>::operator()(vector<T>& aX, const matrix_base<T>& aA, const vector<T>& aB)
{
  // error tolerance
  double error = 1;
  
  // sum of the coefficients and old
  // computed x values
  double csum;
  
  // new value x
  double newx = 0, oldx = 0;
  
  aX.clear();
  aX.reserve(aA.columns(), true);
  
  // continue computing new x values until
  // we are close enough to a solution
  while (error >= m_error_tol)
  {
    // loop through every row of the matrix
    for (unsigned int i = 0; i < aA.rows(); i++)
    {
      csum = 0;
      
      // compute sum of the coefficients time the old
      // x values
      for (unsigned int j = 0; j < aA.columns(); j++)
      {
        if (j != i)
        {
          csum += aA(i, j)*aX[j];
        }
      }
      
      // compute the new x by subtracting the sum from the
      // right-hand side of the equation and dividing by the
      // coefficient corresponding to the current x
      //
      // xi = (an - (a1 + a2 + ... + an-1 (not including ai))) / ai
      newx = (aB[i] - csum) / aA(i, i);
      
      //update computed x value
      if (i == 0)
        oldx = aX[i];
      aX[i] = newx;
    }
    
    // error is the difference in the old and new values
    error = fabs(aX[0] - oldx);
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

