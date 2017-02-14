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

/*
*residualForDislocation
*evaluate residual term at element level for dislocations representing by force dipole
*edge, srew and dislocation loop should be properly described in "IGA_dislocation::mark_plane" 
*/
template <class T, int dim>
void model_dislocation<T, dim>::residualForDislocation(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R, double b)
{
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

template <class T, int dim>
void model_dislocation<T, dim>::residualForPointDefect(NURBSMesh<dim>* _mesh,knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R, dealii::Point<dim> _quadPoints, dealii::Point<dim> strength)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
 	if (std::strcmp(this->bcType,"dislocation")!=0) {printf("Imcompatible Model and BC "); exit(-1);}
  IGAValues<dim> fe_values_temp(_mesh, dim, 0);
   std::vector<std::vector<double> > quadPoints(1);
  quadPoints[0].push_back(_quadPoints[0]); quadPoints[0].push_back(_quadPoints[1]); quadPoints[0].push_back(_quadPoints[2]); quadPoints[0].push_back(2.0);
  std::cout<<"quadPoints"<<quadPoints[0][0]<<quadPoints[0][3]<<std::endl;
  fe_values_temp.reinit(cell, &quadPoints);
  for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
		const unsigned int ck = fe_values.system_to_component_index(dof);
    for (unsigned int q=0; q<1; ++q){
			if (ck==2){ R[dof] += -fe_values_temp.shape_grad(dof, q)[2]*strength[2];}
			else if(ck==1) { R[dof] += -fe_values_temp.shape_grad(dof, q)[1]*strength[1];}
			else if(ck==0) { R[dof] += -fe_values_temp.shape_grad(dof, q)[0]*strength[0];}
		}
	}
}

template class model_dislocation<Sacado::Fad::DFad<double>,1>;
template class model_dislocation<Sacado::Fad::DFad<double>,2>;
template class model_dislocation<Sacado::Fad::DFad<double>,3>;
