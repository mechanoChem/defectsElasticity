/*
*essential functions of model_dislocation
*Users may not change anything here for well-posed initBoundValueProb
*/
#include "model/model_dislocation.h"

template <class T, int dim>
model_dislocation<T, dim>::model_dislocation():model<T,dim>()
{
	std::cout<<"dislocation model generated"<<std::endl;
}


template <class T, int dim>
model_dislocation<T, dim>::~model_dislocation (){}

template <class T, int dim>
void model_dislocation<T, dim>::getResidual(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R,unsigned int currentIteration)
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
	residualForDislocation(cell, fe_values, ULocal, R);
	  if (this->params->getBool("enforceWeakBC")) residualForHighOrderBC(cell,fe_values, ULocal, R);
}

/*
*residualForDislocation
*evaluate residual term at element level for dislocations representing by force dipole
*edge, srew and dislocation loop should be properly described in "IGA_dislocation::mark_plane" 
*/
template <class T, int dim>
void model_dislocation<T, dim>::residualForDislocation(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R)
{
	b=this->params->getDouble("b");
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	  //apply edgedislocation.
	if (std::strcmp(this->bcType,"dislocation")!=0) {printf("Imcompatible Model and BC "); exit(-1);}
  if (std::strcmp(this->bcType,"dislocation")==0){
		//edge dislocation
		if(cell.defectFlags[4]==1){
			for(unsigned int planeID=0; planeID<dim*2; planeID++){
				if(cell.planeFlags[planeID]>0){
					double upperDown;
					if(cell.planeFlags[planeID]==2) upperDown=1;
					else if(cell.planeFlags[planeID]==1) upperDown=-1;
					else {printf("unreadable half plane mark"); exit(-1);};
					
				  for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
				    // double cof=cell.controlPoint_cross_surface_cof[(int)dof/3];// dangerous to use, need to modify IGAbase correspondingly. 
				    double cof=1;
				    for(unsigned int i=0;i<dim;i++){
							if(cell.defectFlags[i]>0){
					  		double mult=cell.defectFlags[i];
					  		const unsigned int ck = fe_values.system_to_component_index(dof) - this->DOF;
					  		for (unsigned int q=0; q<fe_values.n_plane_quadrature_points; ++q){
					   			if( ck==i){
					      		R[dof]+=-fe_values.shape_grad_plane(dof, q, planeID)[ck]*(upperDown*3*b/2)*fe_values.JxW_plane(q, planeID)*cof*mult;
					    		}
					    		else{
					     			R[dof] += -fe_values.shape_grad_plane(dof, q, planeID)[ck]*(upperDown*b/2)*fe_values.JxW_plane(q, planeID)*cof*mult;
					    		}
					  		}
							}
				    }
				  }
				}
			}
		}
		
		if(cell.defectFlags[5]==1){
			for(unsigned int planeID=0; planeID<dim*2; planeID++){
				if(cell.planeFlags[planeID]>0){
					double upperDown;
					if(cell.planeFlags[planeID]==2) upperDown=1;
					else if(cell.planeFlags[planeID]==1) upperDown=-1;
					else {printf("unreadable half plane mark"); exit(-1);};

				  for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
				    double cof=1;
				    for(unsigned int i=0;i<dim;i++){
							if(cell.defectFlags[i]>0){
					  		double mult=cell.defectFlags[i];
					  		const unsigned int ck = fe_values.system_to_component_index(dof) - this->DOF;
					  		for (unsigned int q=0; q<fe_values.n_plane_quadrature_points; ++q){
					   			if( ck==planeID/2){
					      		R[dof] += -fe_values.shape_grad_plane(dof, q, planeID)[i]*(upperDown*b/2)*fe_values.JxW_plane(q, planeID)*cof*mult;
					    		}
					   			else if( ck==i){
					      		R[dof] += -fe_values.shape_grad_plane(dof, q, planeID)[planeID/2]*(upperDown*b/2)*fe_values.JxW_plane(q, planeID)*cof*mult;
					    		}
					  		}
							}
				    }
				  }
				}
			}
		}		
	}
}

template class model_dislocation<Sacado::Fad::DFad<double>,1>;
template class model_dislocation<Sacado::Fad::DFad<double>,2>;
template class model_dislocation<Sacado::Fad::DFad<double>,3>;