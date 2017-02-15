/*
*IGA_structure constructor and destructor
*/

#include "initBoundValueProb/IGA_structure.h"

using namespace std;
template <int dim>
IGA_structure<dim>::IGA_structure (NURBSMesh<dim>& _mesh, parametersClass<dim>& _params):IGA<dim>(_mesh,_params)
{
	IGA<dim>::numIncrements=1; IGA<dim>::currentIncrement=0;
}
template <int dim>
IGA_structure<dim>::~IGA_structure(){}

template class IGA_structure<1>;
template class IGA_structure<2>;
template class IGA_structure<3>;
