/*
*model_structure<T, dim>::run
*overwrite model<T, dim>::run to include residualForNewmmanBC
*User may use it directly
*/

#include "model/model_structure.h"
template <class T, int dim>
void model_structure<T, dim>::getResidual(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R,unsigned int currentIteration)
{
	//clean R check 

	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
    if (abs(R[dof].val())>1.0e-13){
    	printf("**************Residual is contaminated**************. Value: %12.4e\n", R[dof].val());
			printf("\n"); exit(-1);
    }
	}
	this->iteration=currentIteration;
	
	residualForMechanics(cell, fe_values, ULocal, R);
	residualForNewmmanBC(cell, fe_values, ULocal, R);
	//residualForBodyforce(cell, fe_values, ULocal, R);
	  if (this->params->getBool("enforceWeakBC")) residualForHighOrderBC(cell,fe_values, ULocal, R);
}

template class model_structure<Sacado::Fad::DFad<double>,1>;
template class model_structure<Sacado::Fad::DFad<double>,2>;
template class model_structure<Sacado::Fad::DFad<double>,3>;