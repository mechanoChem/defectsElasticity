/*
*IGA_dislocation<dim>::assemble_system_interval
*this function assembles jacobian and righ hand side 
*use "dislocationModel" defined in "run"
*User may use it directly
*/

#include "initBoundValueProb/IGA_dislocation.h"

using namespace std;

template <int dim>
void IGA_dislocation<dim>::assemble_system_interval (const typename std::vector<knotSpan<dim> >::iterator &begin, const typename std::vector<knotSpan<dim> >::iterator &end){
  //element loop
  IGAValues<dim> fe_values_base(IGA<dim>::mesh, dim, 2);
  for (typename std::vector<knotSpan<dim> >::iterator cell=begin; cell<end; cell++){
    fe_values_base.reinit(*cell);
    IGAValues<dim>& fe_values=fe_values_base;
    //IGAValues<dim>* fe_values=cellValues[cell->id];
    unsigned int n_q_points= fe_values.n_quadrature_points;
    unsigned int dofs_per_cell=fe_values.dofs_per_cell;
    denseMatrix local_matrix(dofs_per_cell, dofs_per_cell);
    denseVector local_rhs(dofs_per_cell);
    //AD variables
    dealii::Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); dealii::Table<1, double > ULocalConv(dofs_per_cell);
    dealii::Table<1, double> ULocalTemp(dofs_per_cell); 
    for (unsigned int i=0; i<dofs_per_cell; ++i){
      ULocal[i]=U(cell->local_dof_indices[i]);
      ULocalTemp[i]=U(cell->local_dof_indices[i]);
      ULocal[i].diff (i, dofs_per_cell);
      ULocalConv[i]= Un(cell->local_dof_indices[i]);
    }
	
    dealii::Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); //R.fill(0.0);
    
    for(unsigned int i=0; i<dofs_per_cell; ++i) R[i]=0.0; 
		//gradient elasticity
		dislocationModel->getResidualIni(R, IGA<dim>::currentIteration);
		dislocationModel->residualForMechanics(*cell, fe_values, ULocal, R);
		
		//make vectors for quadPoints and strength for doing multiple pointDefects
	  if (cell->defectFlags[3]==1){
			dealii::Point<dim> quadPoints=IGA<dim>::params.getPoint("DefectQuad1");
			dealii::Point<dim> strength=IGA<dim>::params.getPoint("DefectStrength1");
			dislocationModel->residualForPointDefect(IGA<dim>::mesh, *cell, fe_values, ULocal, R, quadPoints,strength);
		}
		//if (this->params->getBool("enforceWeakBC")) residualForHighOrderBC(cell,fe_values, ULocal, R);
		
		
		
    //Residual(R) and Jacobian(R')
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      for (unsigned int j=0; j<dofs_per_cell; ++j){
				// R' by AD
				local_matrix(i,j)= R[i].fastAccessDx(j);
      }
      local_rhs(i) = -R[i].val();
    }
	
    //Global Assembly
    IGA<dim>::assembler_lock.acquire();
    for (unsigned int i=0; i<dofs_per_cell; ++i){
      for (unsigned int j=0; j<dofs_per_cell; ++j){
				IGA<dim>::system_matrix(cell->local_dof_indices[i], cell->local_dof_indices[j])+=local_matrix(i,j);
      }
      IGA<dim>::system_rhs(cell->local_dof_indices[i]) += local_rhs(i);
    }
    IGA<dim>::assembler_lock.release();
  }
}

template class IGA_dislocation<1>;
template class IGA_dislocation<2>;
template class IGA_dislocation<3>;
