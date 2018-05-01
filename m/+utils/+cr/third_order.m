function res=third_order(dvvv,dvv,dv,vzzz,vzz,vz,options)
% third_order -- third-order multivariate chain rule
%
% ::
%
%
%   res=third_order(dvvv,dvv,dv,vzzz,vzz,vz)
%
%   res=third_order(dvvv,dvv,dv,vzzz,vzz,vz,options)
%
% Args:
%
%    - **dvvv** [nd x nv^3 matrix]: matrix of third derivatives of the d
%    function with respect to its locations. The derivatives are unfolded
%    columnwise
%
%    - **dvv** [nd x nv^2 matrix]: matrix of second derivatives of the d
%    function with respect to its locations. The derivatives are unfolded
%    columnwise
%
%    - **dv** [nd x nv matrix]: jacobian of function with respect to the
%    locations of its arguments
%
%    - **vzzz** [nv x nz^3 matrix]: third derivatives of the locations with
%    respect to the variables to differentiate. The derivatives are unfolded
%    columnwise
%
%    - **vzz** [nv x nz^2 matrix]: second derivatives (hessian) of the
%    locations with respect to the variables to differentiate. The derivatives
%    are unfolded columnwise
%
%    - **vz** [nv x nz matrix]: jacobian of the locations with respect to the
%    variables to differentiate
%
%    - **options** [empty|struct]: When not empty, options is a structure with
%    fields:
%      - **large** [true|{false}] if true, a computation explicitly using the
%      kronecker product is avoided.
%      - **multiply** [true|{false}]: if true, explicit omega matrices are
%      constructed and then multiplied to other matrices to sum the
%      permutations. Else, a functional form is used instead.
%
% Returns:
%    :
%
%    - **res** [nd x nz^3]: output matrix
%
% Note:
%
% Example:
%
%    See also:

if nargin<7
    options=[];
end

default_options={
    'large',false,@(x)islogical(x),'large must be a logical'
    'multiply',false,@(x)islogical(x),'multiply must be a logical'
    };

if isempty(options)
    options=cell2struct(default_options(:,2),default_options(:,1),1);
else
    if ~isstruct(options)
        error('options must be a structure or empty')
    end
    options=parse_arguments(default_options,options);
end

nz=size(vz,2);
nd=size(dv,1);

template=sparse(nd,nz^3);
res=template;

is_computable=@utils.cr.is_computable;

is_dv=is_computable(dv);
is_dvv=false;
is_dvvv=false;

if is_dv
    is_dvv=is_computable(dvv);
    if is_dvv
        is_dvvv=is_computable(dvvv);
    end
    if is_computable(vz)
            res=res+dvvv_vz_vz_vz();
        if is_computable(vzz)
            res=res+dvv_vz_vzz();
            if is_computable(vzzz)
                res=res+dv*vzzz;
            end
        end
    end
end

    function res=dvvv_vz_vz_vz()
        res=template;
        if is_dvvv
            if options.large
                res=utils.kronecker.A_times_kron_Q1_Qk(...
                    dvvv,vz,vz,vz);
            else
                res=dvvv*utils.kronecker.kronall(vz,vz,vz);
            end
        end
    end

    function res=dvv_vz_vzz()
        res=template;
        if is_dvv
            if options.large
                tmp=utils.kronecker.A_times_kron_Q1_Qk(dvv,vz,vzz);
            else
                tmp=dvv*kron(vz,vzz);
            end
            if options.multiply
                res=tmp*utils.cr.omega(nz,1);
            else
                res=utils.cr.dv_vz_omega(tmp,nz,1);
            end
        end
    end

end