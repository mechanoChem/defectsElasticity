template <class T, int dim>
void model_structure<T, dim>::residualForBodyforce(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R)
{

}

template class model_structure<Sacado::Fad::DFad<double>,1>;
template class model_structure<Sacado::Fad::DFad<double>,2>;
template class model_structure<Sacado::Fad::DFad<double>,3>;
