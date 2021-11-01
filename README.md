# weaving-geodesic-foliations
Research Code for Processing Vector Fields on Manifolds

This is the "relax-field" app developed as part of the weaving geodesic foliations paper. https://dl.acm.org/doi/10.1145/3306346.3323043

The app interface is a bit unintuitive to use at the moment, will try to document the intricacies a bit at some point.  

The high level is that this code implements a solver for finding geodesic vector fields on branched covers of 2 manifolds, 
and then additionally includes code for integrating these geodesic fields into geodesic foliations using a nonlinear gauss-newton optimization, and then extracts level sets to generate weaving patterns which serve as inputs for the forward elastic rod simulator implemented here: https://github.com/evouga/RibbonSim
