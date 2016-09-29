#ifndef bFitting_H_
#define bFitting_H_
#include "igaClasses.h"
#include "base/table.h"
#include "solutionClasses.h"
#include "parameters.h"

template <int dim>
double displacementFitting(NURBSMesh<dim>* mesh,denseVector U,parametersClass* params){	
	const double pi=3.1415926;
  double lambda= params->getDouble("lambda");
  double mu=     params->getDouble("mu");
	double Po=lambda/(2*(lambda+mu));
	double R1=0;
	double R2=0;
	double b;
	IGAValues<dim> fe_values_base(mesh, dim, 2);

	std::cout<<"disPlaementFitting"<<std::endl;

	for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){

    fe_values_base.reinit(*cell);
    IGAValues<dim> fe_values=fe_values_base;

    unsigned int n_q_points= fe_values.n_quadrature_points;
    unsigned int dofs_per_cell=fe_values.dofs_per_cell;
		unsigned int controlpoints_per_cell=cell->controlPoints.size();
		//dealii::Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell);
		dealii::Table<1,double> ULocal(dofs_per_cell);
		dealii::Table<2,double> disp_num(n_q_points, dim);
		dealii::Table<2,double> disp_anl_Db(n_q_points, dim);
	
    for (unsigned int i=0; i<dofs_per_cell; ++i){
      ULocal[i]=U(cell->local_dof_indices[i]);
    }

		//numerical displacements
		for (unsigned int q=0; q<n_q_points; ++q){
			for (unsigned int i=0;i<dim; ++i) {
				disp_num[q][i]=0.0;
			}		
	    for (unsigned int dof=0; dof<dofs_per_cell; ++dof){
	      const unsigned int ck = fe_values.system_to_component_index(dof);
				disp_num[q][ck]+=fe_values.shape_value(dof, q)*ULocal[dof];
	    }
		}
		
		//analytial Du,b
		for(unsigned int q=0;q<n_q_points;++q){
			//physical coords
      double x=0;
			double y=0;
			double z=0;
      for (unsigned int basisIndex=0; basisIndex<controlpoints_per_cell; basisIndex++) {
					x+=fe_values.shape_value(basisIndex*dim, q)*cell->controlPoints[basisIndex]->coords[0];
					y+=fe_values.shape_value(basisIndex*dim+1, q)*cell->controlPoints[basisIndex]->coords[1];
					z+=fe_values.shape_value(basisIndex*dim+2, q)*cell->controlPoints[basisIndex]->coords[2];				
			}
		x=x-0.5;
		y=y-0.5;
		z=z-0.5;
			disp_anl_Db[q][0]=1/(2*pi)*(std::atan(y/x)+x*y/(2*(1-Po)*(std::pow(x,2)+std::pow(y,2))));
			disp_anl_Db[q][1]=-1/(2*pi)*((1-2*Po)/4/(1-Po)*std::log(std::pow(x,2)+std::pow(y,2))+(std::pow(x,2)-std::pow(y,2)/4/(1-Po)/(std::pow(x,2)+std::pow(y,2))));
			disp_anl_Db[q][2]=0;
		}
		//assembling

		for(unsigned int q=0;q<n_q_points;++q){
			for(unsigned int i=0;i<dim;i++){
				R1+=disp_num[q][ i]*disp_anl_Db[q][i]*fe_values.JxW(q);
				R2+=disp_anl_Db[q][i]*disp_anl_Db[q][i]*fe_values.JxW(q);
			}
		}
	
	}
	b=R1/R2;
	std::cout<<"R1="<<R1<<"  R2="<<R2<<"  b="<<b<<std::endl;
	return b;

}

template <int dim>
void Fill_analytical_displcment(solutionClass<dim>& displacement_anl,NURBSMesh<dim>* mesh,parametersClass* params){
	const double pi=3.1415926;
	double b=			 params->getDouble("b");
  double lambda= params->getDouble("lambda");
  double mu=     params->getDouble("mu");
	double Po=lambda/(2*(lambda+mu));
	
  //fill displacement Vector
    for (typename std::vector<controlPoint<dim> >::iterator controlpoint=mesh->controlPointVector.begin(); controlpoint<mesh->controlPointVector.end(); controlpoint++){
			std::vector<double> coords(controlpoint->coords);
      			double x=coords[0]-0.5;
			double y=coords[1]-0.5;
			double z=coords[2]-0.5;
			
			displacement_anl(controlpoint->id,0)=b/(2*pi)*(std::atan(y/x)+x*y/(2*(1-Po)*(std::pow(x,2)+std::pow(y,2))));
			displacement_anl(controlpoint->id,1)=-b/(2*pi)*((1-2*Po)/4/(1-Po)*std::log(std::pow(x,2)+std::pow(y,2))+(std::pow(x,2)-std::pow(y,2)/4/(1-Po)/(std::pow(x,2)+std::pow(y,2))));
			displacement_anl(controlpoint->id,2)=0.0;
  }
}

template <int dim>
double get_L2_norm(NURBSMesh<dim>* mesh,denseVector U,parametersClass* params){
	const double pi=3.1415926;
  double lambda= params->getDouble("lambda");
  double mu=     params->getDouble("mu");
	double Po=lambda/(2*(lambda+mu));
	double b=			 params->getDouble("b");
	IGAValues<dim> fe_values_base(mesh, dim, 2);
	double L2_norm;
	double norm1=0;
	double norm2=0;

	std::cout<<"calculating L2_norm"<<std::endl;

	for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){

    fe_values_base.reinit(*cell);
    IGAValues<dim> fe_values=fe_values_base;

    unsigned int n_q_points= fe_values.n_quadrature_points;
    unsigned int dofs_per_cell=fe_values.dofs_per_cell;
		unsigned int controlpoints_per_cell=cell->controlPoints.size();
		//dealii::Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell);
		dealii::Table<1,double> ULocal(dofs_per_cell);
		dealii::Table<2,double> disp_num(n_q_points, dim);
		dealii::Table<2,double> disp_anl(n_q_points, dim);
	
    for (unsigned int i=0; i<dofs_per_cell; ++i){
      ULocal[i]=U(cell->local_dof_indices[i]);
    }

		//numerical displacements
		for (unsigned int q=0; q<n_q_points; ++q){
			for (unsigned int i=0;i<dim; ++i) {
				disp_num[q][i]=0.0;
			}		
	    for (unsigned int dof=0; dof<dofs_per_cell; ++dof){
	      const unsigned int ck = fe_values.system_to_component_index(dof);
				disp_num[q][ck]+=fe_values.shape_value(dof, q)*ULocal[dof];
	    }
		}
		
		//analytial Du,b
		for(unsigned int q=0;q<n_q_points;++q){
			//physical coords
      double x=0;
			double y=0;
			double z=0;
      for (unsigned int basisIndex=0; basisIndex<controlpoints_per_cell; basisIndex++) {
					x+=fe_values.shape_value(basisIndex*dim, q)*cell->controlPoints[basisIndex]->coords[0];
					y+=fe_values.shape_value(basisIndex*dim+1, q)*cell->controlPoints[basisIndex]->coords[1];
					z+=fe_values.shape_value(basisIndex*dim+2, q)*cell->controlPoints[basisIndex]->coords[2];				
			}
		x=x-0.5;
		y=y-0.5;
		z=z-0.5;
			disp_anl[q][0]=b/(2*pi)*(std::atan(y/x)+x*y/(2*(1-Po)*(std::pow(x,2)+std::pow(y,2))));
			disp_anl[q][1]=-b/(2*pi)*((1-2*Po)/4/(1-Po)*std::log(std::pow(x,2)+std::pow(y,2))+(std::pow(x,2)-std::pow(y,2)/4/(1-Po)/(std::pow(x,2)+std::pow(y,2))));
			disp_anl[q][2]=0;
		}
		//assembling

		for(unsigned int q=0;q<n_q_points;++q){
			double tem1=0;
			double tem2=0;
			for(unsigned int i=0;i<dim;i++){
				tem1+=std::pow(disp_anl[q][i]-disp_num[q][i],2);
			  tem2+=std::pow(disp_num[q][i],2);
			}
			norm1+=tem1*fe_values.JxW(q);
			norm2+=tem2*fe_values.JxW(q);
		}
	
	}
	L2_norm=norm1/norm2;
	std::cout<<"norm1="<<norm1<<"  norm2="<<norm2<<"  L_norm="<<L2_norm<<std::endl;
	return L2_norm;	
}




#endif 
