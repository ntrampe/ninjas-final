//
//  Filename:     pdeBounds.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                pde class
//


template <class T>
pdeBounds<T>::pdeBounds()
{
  
}


template <class T>
pdeBounds<T>::pdeBounds(T (*aXLower)(T), T (*aXUpper)(T), T (*aYLower)(T), T (*aYUpper)(T), point<T> aBounds)
{
  
}


template <class T>
pdeBounds<T>::~pdeBounds()
{
  
}


template <class T>
T pdeBounds<T>::operator()(const T aX, const T aY)
{
  if (aX == m_bounds.x())
  {
    return m_xLower(aY);
  }
}

