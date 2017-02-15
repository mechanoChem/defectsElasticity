/*
//Mark special face flags
//We have 7 slots in planeFlags. 0-5 denotes 6 planes. 6: denote plane inside cell or just cell surface. 
//planeFlags[i]=1:i:normal vector of plane;  1->negtive half plane; =2->postive plane;i=0-5
//planeFlags[6]=1 mark plane inside of cell. planeFlags[6]=0 (default) mark plane on cell surface.
//using Sp_planeQuad_Point to specify the position of plane inside of cell.
// 
//We have 6 slots in defectFlags. 0-2 for direction of burger vector. 3:PointDefects 4:edge, 5:screw
//defectFlags[i]=n: i:direction of burger vector, n:number of half planes overlaped. i<=2;
//
// 
*/

#include "initBoundValueProb/IGA_dislocation.h"

using namespace std;

template <int dim>
void IGA_dislocation<dim>::mark_plane(){
  for (typename std::vector<knotSpan<dim> >::iterator cell=IGA<dim>::mesh->knotSpanVector.begin(); cell<IGA<dim>::mesh->knotSpanVector.end(); cell++){
    
		//edge dislocation at center
		if(cell->endKnots[0][1]>0.5 and cell->endKnots[0][0]<0.5 and cell->endKnots[1][0]>=0.5) {
			cell->planeFlags[0*2+1]=2; cell->planeFlags[6]=1; cell->Sp_planeQuad_Point=0;
			cell->defectFlags[4]=1;
		  cell->defectFlags[0]=cell->defectFlags[1]+1;   
		}
		if(cell->endKnots[0][1]>0.5 and cell->endKnots[0][0]<0.5 and cell->endKnots[1][0]<=0.5) {
			cell->planeFlags[0*2+1]=1; cell->planeFlags[6]=1; cell->Sp_planeQuad_Point=0;
			cell->defectFlags[4]=1;
		  cell->defectFlags[0]=cell->defectFlags[1]+1;   
		}
  }
}

template class IGA_dislocation<1>;
template class IGA_dislocation<2>;
template class IGA_dislocation<3>;
