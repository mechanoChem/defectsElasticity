/*
*run
*entrance of run IGA_structure
*User may use it directly
*/
#include "initBoundValueProb/IGA_structure.h"

using namespace std;

template <int dim>
void IGA_structure<dim>::run (){
  model_structure<Sacado::Fad::DFad<double>, dim> _model_structure;
	structureModel=&_model_structure;
	structureModel->reinit(IGA<dim>::params);
	
  IGA<dim>::setup();
  IGA<dim>::mark_boundaries();//normally just using default
	double values[]={0};
  for (unsigned int i=0; i<sizeof(values)/sizeof(double); ++i){
    IGA<dim>::params.setDouble("l",values[i]);
    IGA<dim>::apply_initial_values();//default(Un=0)
    for (IGA<dim>::currentIncrement=1; IGA<dim>::currentIncrement<=IGA<dim>::numIncrements; ++IGA<dim>::currentIncrement){
		        apply_boundary_conditions();
			IGA<dim>::solve(); 
    }
	}
    output(0);
}

template class IGA_structure<1>;
template class IGA_structure<2>;
template class IGA_structure<3>;
