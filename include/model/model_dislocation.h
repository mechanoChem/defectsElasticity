/*
*base class for dislocation modeling
*derived from model
*/

#ifndef dislocationClass_h
#define dislocationClass_h

#include "model/model.h"

template <class T, int dim>
class model_dislocation : public model<T, dim>
{
 public:
 	/*
 	*model_dislocation constructor and destructor
 	*/
  model_dislocation(); // Class constructor 
  ~model_dislocation(); //Class destructor
	
	/*
	*evaluate residual term at element level for dislocations representing by force dipole
	*edge, srew and dislocation loop should be properly described in "IGA_dislocation::mark_plane" 
	*/
	void residualForDislocation(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R, double b);
	void residualForPointDefect(NURBSMesh<dim>* _mesh, knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R, dealii::Point<dim> _quadPoints, dealii::Point<dim> strength);
	
	/*
	*burger vector of dislcoation
	*/
	double b;


};

#endif
