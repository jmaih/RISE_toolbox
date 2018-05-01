function [siz,t,z,parts,order_var,inv_order_var]=solution_topology(...
lead_lag_incidence,...
exo_nbr,... number of shocks
kfuture) % number of shocks periods beyond the current
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

    if nargin<3
        kfuture=0;
    end

[parts]=utils.solve.partition_variables(lead_lag_incidence);

static=parts.static;
pred=parts.pred;
both=parts.both;
frwrd=parts.frwrd;
order_var=parts.order_var;
inv_order_var=parts.inv_order_var;

siz=struct('ns',numel(static),...
    'np',numel(pred),...
    'nb',numel(both),...
    'nf',numel(frwrd),...
    'ne',exo_nbr);

% z-locations
%------------
zp=1:siz.np;
zb=siz.np+(1:siz.nb);
z=struct('p',zp,...
    'b',zb,...
    'sig',siz.np+siz.nb+1,...
    'e_0',siz.np+siz.nb+1+(1:siz.ne),...
    'e_plus',siz.np+siz.nb+1+siz.ne+(1:kfuture*siz.ne),...
    'pb',[zp,zb]);

% t-locations
%------------
tp=siz.ns+(1:siz.np);
tb=siz.ns+siz.np+(1:siz.nb);
tf=siz.ns+siz.np+siz.nb+(1:siz.nf);
t=struct('s',1:siz.ns,...
    'p',tp,...
    'b',tb,...
    'f',tf,...
    'pb',[tp,tb],...
    'bf',[tb,tf]);

end
