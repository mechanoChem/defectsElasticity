/*
*IGA_structure<dim>::apply_boundary_conditions
*this function generates apply dirchlet boundary condition
*
*model_structure<T, dim>::residualForNeummanBC
*apply neumman boundary condition at element level
*
*User may change it accordingly 
*/

#include "initBoundValueProb/IGA_structure.h"

using namespace std;

template <int dim>
void IGA_structure<dim>::apply_boundary_conditions(){
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

template class IGA_structure<1>;
template class IGA_structure<2>;
template class IGA_structure<3>;


#include "model/model_structure.h"
//comment
/*This function does this
* and this
* and this....
*/
template <class T, int dim>
void model_structure<T, dim>::residualForNeummanBC(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R)
{
	double load= model<T, dim>::params->getDouble("load"); 
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if (cell.boundaryFlags[faceID]>0){      
      //loop over boundary faces
      for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
				const unsigned int ck = fe_values.system_to_component_index(dof) - model<T, dim>::DOF;
				for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
	  			if (std::strcmp(model<T, dim>::bcType,"tension")==0){
	    			//tension along dim-1 direction
	    			if ((cell.boundaryFlags[faceID]==4) and (ck==1)){
	      			R[dof] += -fe_values.shape_value_face(dof, q, faceID)*load*fe_values.JxW_face(q, faceID);
	    			}
	  			}
	  		}
			} 
    } 
  }
}

template class model_structure<Sacado::Fad::DFad<double>,1>;
template class model_structure<Sacado::Fad::DFad<double>,2>;
template class model_structure<Sacado::Fad::DFad<double>,3>;
