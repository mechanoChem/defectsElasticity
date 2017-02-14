/*
*output
*this function generates vtk file
*User may use it directly
*/

#include "initBoundValueProb/IGA_dislocation.h"

using namespace std;

template <int dim>
void IGA_dislocation<dim>::output (unsigned int _cycle){
	solutionClass<dim> displacement(*IGA<dim>::mesh, NODAL, VECTOR, std::string("u"));
	//fill output vector
  for (typename std::vector<controlPoint<dim> >::iterator i=IGA<dim>::mesh->controlPointVector.begin(); i<IGA<dim>::mesh->controlPointVector.end(); i++){
    for (unsigned int a=0; a<dim; a++){
			displacement(i->id, a)=IGA<dim>::Un(i->id*IGA<dim>::mesh->dofPerControlPoint + a);
    }
  }
	IGA<dim>::outputVariables.push_back(&displacement);
  char fileName[200];
  std::sprintf (fileName, "%s/output",this->params.getString("output_path").c_str());
  std::cout<<"begin write mesh"<<std::endl;
  std::vector<int> outputGridSize(dim,10);
  writeMesh<dim>(fileName, _cycle, IGA<dim>::mesh, outputGridSize, IGA<dim>::outputVariables);
}

template class IGA_dislocation<1>;
template class IGA_dislocation<2>;
template class IGA_dislocation<3>;
