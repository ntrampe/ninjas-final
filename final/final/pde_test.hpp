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
  return sin(M_PI/2.0 - 5*aY) + sin(0.5*aY) + 2*sin(2*aY);
}


template <class T>
T pde_test<T>::xUpper(T aY) const
{
  return sin(M_PI/2.0 - 5*aY) + sin(0.5*aY) + 2*sin(2*aY) - 3;
}


template <class T>
T pde_test<T>::yLower(T aX) const
{
  return sin(M_PI/2.0 - 5*aX) + sin(0.5*aX) + 2*sin(2*aX);
}


template <class T>
T pde_test<T>::yUpper(T aX) const
{
  return sin(M_PI/2.0 - 5*aX) + sin(0.5*aX) + 2*sin(2*aX) - 3;
}

