//
//  Filename:     gauss_elim.hpp
//  Programmer:   James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - A Parameterized Gaussian Elimination Class
//
//  Description:  Here are the definitions of the functions declared in the gauss_elim class
//
//


template <class T>
bool gauss_elim<T>::operator()(vector<T>& aX, const matrix_base<T>& aA, const vector<T>& aB)
{
	matrix<T> augMat(aA);

	vector<T> b(aB);

	aX.clear();
	aX.reserve(augMat.rows(), true);

	int rows = static_cast<int>(augMat.rows());
	int columns = static_cast<int>(augMat.columns());
	int j, col, row, max_row, diag;
	T max, tmp;

	for (diag = 0; diag < columns; diag++)
	{
		max_row = diag, max = augMat(diag, diag);

		for (row = diag + 1; row < rows; row++)
			if ((tmp = fabs(augMat(row, diag))) > max)
				max_row = row, max = tmp;

		augMat.switchRows(diag, max_row);
		T tmp = b[diag];
		b[diag] = b[max_row];
		b[max_row] = tmp;

		for (row = diag + 1; row < rows; row++)
		{
			tmp = augMat(row, diag) / augMat(diag, diag);
			for (col = diag+1; col < columns; col++)
				augMat(row, col) -= tmp * augMat(diag, col);
			augMat(row, diag) = 0;
			b[row] -= tmp * b[diag];
		}
	}

	for (row = columns - 1; row >= 0; row--)
	{
		tmp = b[row];
		for (j = rows - 1; j > row; j--)
			tmp -= aX[j] * augMat(row, j);
		aX[row] = tmp / augMat(row, row);
	}

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

