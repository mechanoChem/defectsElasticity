#include "initBoundValueProb/IGA_dislocation.h"

using namespace std;

template <int dim>
void IGA_dislocation<dim>::run (){
  model_dislocation<Sacado::Fad::DFad<double>, dim> _model_dislocation;
	dislocationModel=&_model_dislocation;
	dislocationModel->reinit(IGA<dim>::params);
	
  IGA<dim>::setup();
  IGA<dim>::mark_boundaries();//normally just using default
	mark_plane();
	double values[]={0};
  for (unsigned int i=0; i<sizeof(values)/sizeof(double); ++i){
    IGA<dim>::params.setDouble("l",values[i]);
    IGA<dim>::apply_initial_values();//default(Un=0)
    //output(0);
    for (IGA<dim>::currentIncrement=1; IGA<dim>::currentIncrement<=IGA<dim>::numIncrements; ++IGA<dim>::currentIncrement){
			IGA<dim>::apply_boundary_conditions();
			IGA<dim>::solve(); 
			//L2_norm=get_L2_norm( mesh, U, params);
			//output(currentIncrement);
    }
	}
	output(0);
}

template class IGA_dislocation<1>;
template class IGA_dislocation<2>;
template class IGA_dislocation<3>;

