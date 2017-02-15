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
template <class T, int dim>
void model_structure<T, dim>::residualForNewmmanBC(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R)
{
	double load= model<T, dim>::params->getDouble("load"); 
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
//example for Newmman boundary condition; no need for defect paper.
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if (cell.boundaryFlags[faceID]>0){      
      //loop over boundary faces
      for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
				const unsigned int ck = fe_values.system_to_component_index(dof) - model<T, dim>::DOF;
				for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
					if (std::strcmp(model<T, dim>::bcType,"torsion")==0){
						if (dim<3) throw "torsion B.C only valid for 3D problem";
						//torque on z=1 face
						if ((cell.boundaryFlags[faceID]==2*dim) and (ck<=dim-2)){
							unsigned int fq=fe_values.quadratureMap[faceID][q]; //quadrature point index for face quadrature points
						  double x=fe_values.quadPointLocations[fq][0]-0.5, y=fe_values.quadPointLocations[fq][1]-0.5;
						  double r=std::sqrt(x*x+y*y);
						  double theta=std::atan2(y,x);
						  double loadX=load*r*std::sin(theta);
						  double loadY=-load*r*std::cos(theta);
						  if (ck==0) {
								R[dof] += -fe_values.shape_value_face(dof, q, faceID)*loadX*fe_values.JxW_face(q, faceID);
						  }
						  else if (ck==1) {
								R[dof] += -fe_values.shape_value_face(dof, q, faceID)*loadY*fe_values.JxW_face(q, faceID);
						  }
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
