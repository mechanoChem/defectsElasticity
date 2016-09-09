defectsElasticity code: structure and defects modeling based on openIGA 

from the Computational Physics Group at the University of Michigan, adviser: Prof. Krishna Garikipati

current lead developer: Zhenlin Wang

List of contributors:

Zhenlin Wang

Shiva Rudraraju


Overview
=======================================================================
defectsElasticity code is an isogeometric analysis based code, for solving the partial differential equations. Currently it includes two initBoundValueProblems.

1.structure problem [S. Rudraraju, A. Van der Ven, K. Garikipati, “Three dimensional iso-geometric solutions to general boundary value problems of Toupin's theory of gradient elasticity at finite strains“, Computer Methods in Applied Mechanics and Engineering, 2014]

2.defects represents by force dipole [Z, Wang, S. Rudraraju, K. Garikipati, “A three dimensional field formulation, and isogeometric solutions to point and line defects using Toupin's theory of gradient elasticity at finite strains”, Journal of the Mechanics and Physics of Solids, 2016]

The code is based on openIGA lib [https://github.com/mechanoChem/openIGA]


Version information
=======================================================================
This is version 0.1, the intial release of the code.


License
=======================================================================
GNU LESSER GENERAL PUBLIC LICENSE. Please see the file LICENSE for details.


Acknowledgements
=======================================================================
This code has been developed under the support of the following:

1. NSF CDI Type I Grant: CHE1027729 “Meta-Codes for Computational Kinetics”

2. NSF DMREF grant: DMR1436154 “DMREF: Integrated Computational Framework for Designing Dynamically Controlled Alloy-Oxide Heterostructures”.


Usage
=======================================================================
1) download and install openIGA [https://github.com/mechanoChem/openIGA]

2) download defectElasticity code

3) modify code as necessary; initBoundValueProblems contains templates for defects problem and structure problem

4) modify CmakeLists.txt 

5) $cmake CmakeLists.txt

6) $make run


Reference
=======================================================================
If you write a paer using results obtained with the help of this code,  please consider citing one or more of the following:

1) Z, Wang, S. Rudraraju, K. Garikipati, “A three dimensional field formulation, and isogeometric solutions to point and line defects using Toupin's theory of gradient elasticity at finite strains”, Journal of the Mechanics and Physics of Solids, Vol. 94: 336-361, 2016, doi:10.1016/j.jmps.2016.03.028 [http://arxiv.org/abs/1405.3293]

2) S. Rudraraju, A. Van der Ven, K. Garikipati, “Three dimensional iso-geometric solutions to general boundary value problems of Toupin's theory of gradient elasticity at finite strains”, Computer Methods in Applied Mechanics and Engineering Vol 278: 705-728, 2014 [http://arxiv.org/abs/1404.0094]


data:08/22/2016
