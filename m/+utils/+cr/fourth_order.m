function res=fourth_order(dvvvv,dvvv,dvv,dv,vzzzz,vzzz,vzz,vz,options)
% fourth_order -- fourth-order multivariate chain rule
%
% ::
%
%
%   res=fourth_order(dvvvv,dvvv,dvv,dv,vzzzz,vzzz,vzz,vz)
%
%   res=fourth_order(dvvvv,dvvv,dvv,dv,vzzzz,vzzz,vzz,vz,options)
%
% Args:
%
%    - **dvvvv** [nd x nv^4 matrix]: matrix of fourth derivatives of the d
%    function with respect to its locations. The derivatives are unfolded
%    columnwise
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
%    - **vzzzz** [nv x nz^4 matrix]: fourth derivatives of the locations with
%    respect to the variables to differentiate. The derivatives are unfolded
%    columnwise
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
%    - **res** [nd x nz^4]: output matrix
%
% Note:
%
% Example:
%
%    See also:

if nargin<9
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

template=sparse(nd,nz^4);
res=template;

is_computable=@utils.cr.is_computable;

is_dv=is_computable(dv);
is_dvv=false;
is_dvvv=false;
is_dvvvv=false;

if is_dv
    is_dvv=is_computable(dvv);
    if is_dvv
        is_dvvv=is_computable(dvvv);
        if is_dvvv
            is_dvvvv=is_computable(dvvvv);
        end
    end
    if is_computable(vz)
        res=res+dvvvv_vz_vz_vz_vz();
        if is_computable(vzz)
            res=res+dvv_vzz_vzz();
            res=res+dvvv_vz_vz_vzz();
            if is_computable(vzzz)
                res=res+dvv_vz_vzzz();
                if is_computable(vzzzz)
                    res=res+dv*vzzzz;
                end
            end
        end
    end
end

    function res=dvvvv_vz_vz_vz_vz()
        res=template;
        if is_dvvvv
            if options.large
                res=res+utils.kronecker.A_times_kron_Q1_Qk(...
                    dvvvv,vz,vz,vz,vz);
            else
                res=res+dvvvv*utils.kronecker.kronall(vz,vz,vz,vz);
            end
        end
    end

    function res=dvvv_vz_vz_vzz()
        res=template;
        if is_dvvv
            if options.large
                tmp=utils.kronecker.A_times_kron_Q1_Qk(...
                    dvvv,vz,vz,vzz);
            else
                tmp=dvvv*utils.kronecker.kronall(vz,vz,vzz);
            end
            if options.multiply
                res=tmp*utils.cr.omega(nz,2);
            else
                res=utils.cr.dv_vz_omega(tmp,nz,2);
            end
        end
    end

    function res=dvv_vz_vzzz()
        res=template;
        if is_dvv
            if options.large
                tmp=utils.kronecker.A_times_kron_Q1_Qk(dvv,vz,vzzz);
            else
                tmp=dvv*kron(vz,vzzz);
            end
            if options.multiply
                res=res+tmp*utils.cr.omega(nz,3);
            else
                res=res+utils.cr.dv_vz_omega(tmp,nz,3);
            end
        end
    end

    function res=dvv_vzz_vzz()
        res=template;
        if is_dvv
            if options.large
                tmp=utils.kronecker.A_times_kron_Q1_Qk(dvv,vzz,vzz);
            else
                tmp=dvv*kron(vzz,vzz);
            end
            if options.multiply
                res=tmp*utils.cr.omega(nz,4);
            else
                res=utils.cr.dv_vz_omega(tmp,nz,4);
            end
        end
    end
end