/*
*base class for structure problem modeling
*derived from model
*/

#ifndef structureClass_h
#define structureClass_h

#include "model/model.h"

template <class T, int dim>
class model_structure : public model<T, dim>
{
 public:
	 /*
   *model_structure constructor and destructor
   */
   model_structure(); // Class constructor 
   ~model_structure(); //Class destructor
	
 	/*
 	*function to call all evaluate residual term
 	*/
	void getResidual(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R,unsigned int currentIteration);
	
	/*
	*apply neumman boundary condition at element level
	*/
	void residualForNeummanBC(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R);
	
	/*
	*evaluate residual of body force term
	*/
	//void residualForBodyforce(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R);

};

#endif
