//
//  Filename:     pde_test.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                pde class
//


template <class T>
T pde_test<T>::xLower(T aY) const
{
  return cos(1.0/5*aY);
}


template <class T>
T pde_test<T>::xUpper(T aY) const
{
  return cos(1.0/5*aY);
}


template <class T>
T pde_test<T>::yLower(T aX) const
{
  return cos(1.0/5*aX);
}


template <class T>
T pde_test<T>::yUpper(T aX) const
{
  return cos(1.0/5*aX);
}

