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
bool gauss_seidel<T>::operator()(vector<T>& aX, const matrix_base<T>& aA, const vector<T>& aB)
{
  for (unsigned i = 0; i < aMatrix.size(); i++)
  {
    if (aMatrix[i].size() <= aMatrix.size())
      throw std::invalid_argument("Matrix must contain enough columns");
  }
  
  // error tolerance
  double error = 1;
  
  // sum of the coefficients and old
  // computed x values
  double csum;
  
  // new value x
  double newx;
  
  aX.clear();
  aX.reserve(aA.columns());
  
  // continue computing new x values until
  // we are close enough to a solution
  while (error >= aErrorTolerance)
  {
    // loop through every row of the matrix
    for (unsigned int i = 0; i < aA.row(); i++)
    {
      csum = 0;
      
      // compute sum of the coefficients time the old
      // x values
      for (unsigned int j = 0; j < aMatrix.size(); j++)
      {
        if (j != i)
        {
          csum += a[i][j]*x[j];
        }
      }
      
      // compute the new x by subtracting the sum from the
      // right-hand side of the equation and dividing by the
      // coefficient corresponding to the current x
      //
      // xi = (an - (a1 + a2 + ... + an-1 (not including ai))) / ai
      newx = ((a[i][aMatrix.size()] - csum) / a[i][i]);
      
      // error is the difference in the old and new values
      error = fabs(x[i] - newx);
      
      //update computed x value
      x[i] = newx;
    }
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

