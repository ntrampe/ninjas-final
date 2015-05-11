//
//  Filename:     pde_base.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                pde class
//


template <class T>
pde_base<T>::pde_base()
{
	m_density = 5;
	this->m_bounds.set(0,1);
}


template <class T>
pde_base<T>::pde_base(const size_t aN, point2d<T> aBounds)
{
	m_density = aN;
	m_bounds = aBounds;
}


template <class T>
pde_base<T>::pde_base(const size_t aN, const T aLowerBound, const T aUpperBound)
{
	m_density = aN;
	m_bounds.set(aLowerBound, aUpperBound);
}


template <class T>
pde_base<T>::~pde_base()
{
	m_points.clear();
}


template <class T>
size_t pde_base<T>::density() const
{
	return m_density;
}


template <class T>
point2d<T> pde_base<T>::bounds() const
{
	return m_bounds;
}


template <class T>
T pde_base<T>::lowerBound() const
{
	return m_bounds.x();
}


template <class T>
T pde_base<T>::upperBound() const
{
	return m_bounds.y();
}


template <class T>
void pde_base<T>::setDensity(const size_t aN)
{
	m_density = aN;
}


template <class T>
void pde_base<T>::addPoint(const T aX, const T aY, const T aZ)
{
	point2d<T> p(aX, aY);
	addPoint(p, aZ);
}


template <class T>
void pde_base<T>::addPoint(const point2d<T>& aPoint, const T aZ)
{
	point3d<T> p(aPoint.x(), aPoint.y(), aZ);
	m_points.push_back(p);
}


template <class T>
void pde_base<T>::clearPoints()
{
	m_points.clear();
}


template <class T>
std::string pde_base<T>::matlabOutput(float aAnimationFactor, const bool aDrawLines) const
{
	T iBound = 0;
	T jBound = 0;
	T lowerBound = m_bounds.x();
	T upperBound = m_bounds.y();
	double inc = fabs((upperBound - lowerBound)) / m_density;
	double tolerance = inc / 2.0;
	std::stringstream res, ssX, ssY, ssZ;

	// changing the previous for loops with while loops prevents
	// round-off error

	iBound = lowerBound;

	while (std::abs(iBound-(upperBound+inc)) >= tolerance)
	{
		jBound = lowerBound;

		while (std::abs(jBound-(upperBound+inc)) >= tolerance)
		{
			if (equivalent(jBound, lowerBound) || equivalent(iBound, lowerBound) || equivalent(jBound, upperBound) || equivalent(iBound, upperBound))
			{
				ssX << jBound << ", ";
				ssY << iBound << ", ";
				ssZ << this->operator()(jBound, iBound) << ", ";
			}

			jBound += inc;
		}

		iBound += inc;
	}

	for (size_t i = 0; i < m_points.size(); i++)
	{
		ssX << m_points[i].x() << (i == m_points.size()-1 ? "" : ", ");
		ssY << m_points[i].y() << (i == m_points.size()-1 ? "" : ", ");
		ssZ << m_points[i].z() << (i == m_points.size()-1 ? "" : ", ");
	}

	res << "clear;" << std::endl;
	res << "X = [" << ssX.str() << "];" << std::endl;
	res << "Y = [" << ssY.str() << "];" << std::endl;
	res << "Z = [" << ssZ.str() << "];" << std::endl;
	res << "tri = delaunay(X,Y);" << std::endl;
	res << "fig = trisurf(tri, X, Y, Z);" << std::endl;

	if (aDrawLines)
	{
		res << "set(fig,'LineWidth',0.01);" << std::endl;
	}
	else
	{
		res << "set(fig,'EdgeColor','None');" << std::endl;
	}

	if (aAnimationFactor != 0)
	{
		res << "axis vis3d;" << std::endl;
		res << "axis manual;" << std::endl;
		res << "while ishandle(fig); camorbit(" << aAnimationFactor << ",0.0); drawnow; end" << std::endl;
	}

	return res.str();
}


template <class T>
std::string pde_base<T>::pointsOutput() const
{
	std::stringstream res;

	for (size_t i = 0; i < m_points.size(); i++)
	{
		const point3d<T>* p = &m_points[i];
		res << "u(" << std::setw(10) << p->x() << "," << std::setw(10) << p->y() << ") = " << std::setw(10) << p->z() << std::endl;
	}

	return res.str();
}


template <class T>
T pde_base<T>::operator()(const T aX, const T aY) const
{
	if (equivalent(aX, m_bounds.x()))
	{
		return xLower(aY);
	}
	else if (equivalent(aX, m_bounds.y()))
	{
		return xUpper(aY);
	}
	else if (equivalent(aY, m_bounds.x()))
	{
		return yLower(aX);
	}
	else if (equivalent(aY, m_bounds.y()))
	{
		return yUpper(aX);
	}

	return T(0);
}

