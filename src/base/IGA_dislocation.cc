/*
*IGA_dislocation constructor and destructor
*/
#include "initBoundValueProb/IGA_dislocation.h"

using namespace std;
template <int dim>
IGA_dislocation<dim>::IGA_dislocation (NURBSMesh<dim>& _mesh, parametersClass<dim>& _params):IGA<dim>(_mesh,_params)
{
	IGA<dim>::numIncrements=1; IGA<dim>::currentIncrement=0;
}

template <int dim>
IGA_dislocation<dim>::~IGA_dislocation(){}

template class IGA_dislocation<1>;
template class IGA_dislocation<2>;
template class IGA_dislocation<3>;
