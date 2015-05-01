//
//  Filename:     gauss_elim.hpp
//  Programmer:   James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   7 - A Parameterized Gaussian Elimination Class
//
//  Description:  Here are the definitions of the functions declared in the gauss_elim class
//
//


template <class T>
bool gauss_elim<T>::operator()(vector<T>& aX, const matrix_base<T>& aA, const vector<T>& aB)
{
  //// Gaussian Elimination With Scaled Partial Pivoting

	// augmented matrix
	matrix<T> augMat(aA);

	// multiplier for elimination
	double mult = 0;

	// sum of the coefficients and old
	// computed x values in backwards
	// substitution
	double csum = 0;

	// row of largest element in current column
	size_t maxRow = 0;

	// add b vector to the end of the matrix
	augMat.insertVectorAtColumn(aB, aA.columns());

	// create solution array

	aX.clear();
	aX.reserve(augMat.columns() - 1);

	for (size_t i = 0; i < augMat.columns() - 1; i++)
	{
		aX.push_back(0);
	}

	// elimination
	// loop through all columns
	for (size_t col = 0; col < augMat.columns() - 1; col++)
	{
		// set the largest row to the current pivot
		maxRow = col;

		// check the elements under the pivot to see if one is larger
		for (size_t row = col + 1; row < augMat.rows(); row++)
		{
			if (fabs(augMat(row, col)) > fabs(augMat(maxRow, col)))
			{
				maxRow = row;
			}
		}

		// ensure the largest element in and below the pivot is
		// located at the pivot by switching the rows
		if (maxRow != col)
		{
			augMat.switchRows(maxRow, col);
		}

		// forward elimination
		// (zero-out) all elements underneath the pivot
		for (size_t row = col + 1; row < augMat.rows(); row++)
		{
			if (augMat(row, col) != 0)
			{
				mult = -augMat(col, col) / augMat(row, col);
				augMat.replaceVectorAtRow(mult * augMat.vectorAtRow(row) + augMat.vectorAtRow(col), row);
			}
		}

		// scale all pivots to 1
		augMat.scaleRow(col, 1.0f / augMat(col,col), col);
	}

	// backward substitution
	// solve for the unknown variables from the bottom up
	for (int row = (int)augMat.columns() - 1 - 1; row >= 0; row--)
	{
		csum = 0;

		// compute sum of the coefficients times the old
		// x values
		for (size_t col = row; col < augMat.columns() - 1; col++)
		{
			csum += augMat(row, col) * aX[col];
		}

		// compute the new x by subtracting the sum from the
		// right-hand side of the equation and dividing by the
		// coefficient corresponding to the current x
		//
		// xi = (an - (a1 + a2 + ... + an-1 (not including ai))) / ai
		//
		// where i = row
		aX[row] = (augMat(row, augMat.columns() - 1) - csum) / augMat(row, row);

		// fix rounding errors (e.g. 0 != 0.000000000001)
		// for display purposes
		if (equivalent(aX[row], 0))
		{
			aX[row] = 0;
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

