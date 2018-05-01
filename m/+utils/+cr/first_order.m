function res=first_order(dv,vz,~)
% first_order -- first-order multivariate chain rule
%
% ::
%
%
%   res=first_order(dv,vz)
%
% Args:
%
%    - **dv** [nd x nv matrix]: jacobian of function with respect to the
%    locations of its arguments
%
%    - **vz** [nv x nz matrix]: jacobian of the locations with respect to the
%    variables to differentiate
%
% Returns:
%    :
%
%    - **res** [nd x nz]: output matrix
%
% Note:
%
% Example:
%
%    See also:

nz=size(vz,2);
nd=size(dv,1);

template=sparse(nd,nz);
res=template;

is_computable=@utils.cr.is_computable;

is_dv=is_computable(dv);

if is_dv
    if is_computable(vz)
        res=dv*vz;
    end
end

end