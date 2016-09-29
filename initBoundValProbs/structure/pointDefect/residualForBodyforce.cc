template <class T, int dim>
void model_structure<T, dim>::residualForBodyforce(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R)
{
  if (std::strcmp(bcType,"pointDefect")==0){
    if (cell.endKnots[0][0]<0.5 and cell.endKnots[0][1]>0.5 and cell.endKnots[1][0]<0.5 and cell.endKnots[1][1]>0.5 and cell.endKnots[2][0]<0.5 and cell.end\
Knots[2][1]>0.5){
      IGAValues<dim> fe_values_temp(mesh, dim, 0);
      std::vector<std::vector<double> > quadPoints(1);
      quadPoints[0].push_back(0); quadPoints[0].push_back(0); quadPoints[0].push_back(0); quadPoints[0].push_back(2.0);
      fe_values_temp.reinit(cell, &quadPoints);
      for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
        const unsigned int ck = fe_values.system_to_component_index(dof) - DOF;
        for (unsigned int q=0; q<1; ++q){
          if (ck==2){ R[dof] += -fe_values_temp.shape_grad(dof, q)[2]*(load);}
          else if(ck==1) { R[dof] += -fe_values_temp.shape_grad(dof, q)[1]*(load);}
          else if(ck==0) { R[dof] += -fe_values_temp.shape_grad(dof, q)[0]*(load);}
        }
      }
    }
  }
}

template class model_structure<Sacado::Fad::DFad<double>,1>;
template class model_structure<Sacado::Fad::DFad<double>,2>;
template class model_structure<Sacado::Fad::DFad<double>,3>;
