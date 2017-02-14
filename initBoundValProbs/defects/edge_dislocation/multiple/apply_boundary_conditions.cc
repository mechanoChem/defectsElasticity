#include "initBoundValueProb/IGA_dislocation.h"

using namespace std;

template <int dim>
void IGA_dislocation<dim>::apply_boundary_conditions(){
  //Dirichlet map
  IGA<dim>::dirichletMap.clear(); 
  unsigned int  controlPointDOF=-dim;
  for (typename std::vector<controlPoint<dim> >::iterator controlpoint=IGA<dim>::mesh->controlPointVector.begin(); controlpoint<IGA<dim>::mesh->controlPointVector.end(); controlpoint++){
    std::vector<double> coords(controlpoint->coords);
    controlPointDOF+=dim;
    if( coords[1]==0.0){
      IGA<dim>::dirichletMap[controlPointDOF+0]=0.0;
      IGA<dim>::dirichletMap[controlPointDOF+1]=0.0;
      IGA<dim>::dirichletMap[controlPointDOF+2]=0.0;
    }

    //Apply values to solution vector
    for (std::map<unsigned int, double>::iterator dof=IGA<dim>::dirichletMap.begin(); dof!=IGA<dim>::dirichletMap.end(); dof++){
      IGA<dim>::U(dof->first)=dof->second;
    }
  }
	std::cout<<"finished BC"<<std::endl;
}

template class IGA_dislocation<1>;
template class IGA_dislocation<2>;
template class IGA_dislocation<3>;