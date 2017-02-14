/**
*base class for structure boundary value problem 
*derived from IGA
*/

#ifndef IGA_structureClass_h
#define IGA_structureClass_h

#include "IGA.h"
//method
#include "model/model_structure.h"
using namespace std;
template <int dim>
class IGA_structure: public IGA<dim>
{
public:
	/**
	*IGA_structure constructor and destructor
	*/
  IGA_structure(NURBSMesh<dim>& _mesh, parametersClass<dim>& _params);
  ~IGA_structure();
	
	/**
	*the comments of following functions can be found at IGA_dislocation.h
	*Please see .cc file directly 
	*/
  void run();
  void apply_boundary_conditions();
  void output (unsigned int _cycle);
  void assemble_system_interval (const typename std::vector<knotSpan<dim> >::iterator &begin, const typename std::vector<knotSpan<dim> >::iterator &end);
  
	/**
	*include model type
	*/
	model_structure<Sacado::Fad::DFad<double>,dim>* structureModel;
};

#endif