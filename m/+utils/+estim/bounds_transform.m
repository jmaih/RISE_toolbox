%  transforms constrained parameters into unconstrained ones for use in root
%  finders or optimization routines that do not handle constraints on
%  parameters.
% 
%  returns function handles
%  - tr: transform parameters into unconstrained units
%  - utr: unstransform parameters into original units
% 
%  ------------------------------- tests ----------------------------------
%  [tr,utr]=utils.estim.bounds_transform(0,inf),abs(utr(tr(11))-11)<1e-13
%  [tr,utr]=utils.estim.bounds_transform(-inf,0),abs(utr(tr(-11))+11)<1e-13
%  [tr,utr]=utils.estim.bounds_transform(10,15),abs(utr(tr(11))-11)<1e-13
%  [tr,utr]=utils.estim.bounds_transform(10,inf),abs(utr(tr(11))-11)<1e-13
%  [tr,utr]=utils.estim.bounds_transform(-inf,15),abs(utr(tr(11))-11)<1e-13
%