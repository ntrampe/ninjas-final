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
  this->m_bounds.set(0,1);
}


template <class T>
pde_base<T>::pde_base(point2d<T> aBounds)
{
  m_bounds = aBounds;
}


template <class T>
pde_base<T>::pde_base(const T aLowerBound, const T aUpperBound)
{
  m_bounds.set(aLowerBound, aUpperBound);
}


template <class T>
pde_base<T>::~pde_base()
{
  m_points.clear();
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
void pde_base<T>::addKnownPoint(const T aX, const T aY)
{
  addPoint(aX, aY, this->operator()(aX, aY));
}


template <class T>
void pde_base<T>::addKnownPoint(const point2d<T>& aPoint)
{
  addPoint(aPoint, this->operator()(aPoint.x(), aPoint.y()));
}


template <class T>
std::string pde_base<T>::matlabOutput(float aAnimationFactor, const bool aDrawLines) const
{
  std::stringstream res, ssX, ssY, ssZ;
  
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

