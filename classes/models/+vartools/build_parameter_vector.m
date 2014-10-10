function x1=build_parameter_vector(vdata,a_post,vcov)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


npar=numel(vdata.orig_estim_names);

x1=nan(npar,1);
x1(vdata.estim_locs)=a_post;
    
[SIG,OMG]=vartools.covariance_decomposition(vcov);

pvals=nan(size(vdata.same));
pvals(vdata.same)=SIG(vdata.p1(vdata.same));
pvals(~vdata.same)=diag(OMG(vdata.p1(~vdata.same),vdata.p2(~vdata.same)));

x1(vdata.not_estimated)=pvals;

end