/*
*model_dislocation constructor and destructor
*/
#include "model/model_structure.h"

template <class T, int dim>
model_structure<T, dim>::model_structure():model<T,dim>()
{
	std::cout<<"structure model generated"<<std::endl;
}


template <class T, int dim>
model_structure<T, dim>::~model_structure (){}


template class model_structure<Sacado::Fad::DFad<double>,1>;
template class model_structure<Sacado::Fad::DFad<double>,2>;
template class model_structure<Sacado::Fad::DFad<double>,3>;