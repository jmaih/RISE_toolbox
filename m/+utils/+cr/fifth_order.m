function res=fifth_order(dvvvvv,dvvvv,dvvv,dvv,dv,vzzzzz,vzzzz,vzzz,vzz,vz,options)
% fifth_order -- fifth-order multivariate chain rule
%
% Syntax
% -------
% ::
%
%   res=fifth_order(dvvvvv,dvvvv,dvvv,dvv,dv,vzzzzz,vzzzz,vzzz,vzz,vz)
%
%   res=fifth_order(dvvvvv,dvvvv,dvvv,dvv,dv,vzzzzz,vzzzz,vzzz,vzz,vz,options)
%
% Inputs
% -------
%
% - **dvvvvv** [nd x nv^5 matrix]: matrix of fifth derivatives of the d
% function with respect to its locations. The derivatives are unfolded
% columnwise
%
% - **dvvvv** [nd x nv^4 matrix]: matrix of fourth derivatives of the d
% function with respect to its locations. The derivatives are unfolded
% columnwise
%
% - **dvvv** [nd x nv^3 matrix]: matrix of third derivatives of the d
% function with respect to its locations. The derivatives are unfolded
% columnwise
%
% - **dvv** [nd x nv^2 matrix]: matrix of second derivatives of the d
% function with respect to its locations. The derivatives are unfolded
% columnwise
%
% - **dv** [nd x nv matrix]: jacobian of function with respect to the
% locations of its arguments
%
% - **vzzzzz** [nv x nz^5 matrix]: fifth derivatives of the locations with
% respect to the variables to differentiate. The derivatives are unfolded
% columnwise
%
% - **vzzzz** [nv x nz^4 matrix]: fourth derivatives of the locations with
% respect to the variables to differentiate. The derivatives are unfolded
% columnwise
%
% - **vzzz** [nv x nz^3 matrix]: third derivatives of the locations with
% respect to the variables to differentiate. The derivatives are unfolded
% columnwise
%
% - **vzz** [nv x nz^2 matrix]: second derivatives (hessian) of the
% locations with respect to the variables to differentiate. The derivatives
% are unfolded columnwise
%
% - **vz** [nv x nz matrix]: jacobian of the locations with respect to the
% variables to differentiate
%
% - **options** [empty|struct]: When not empty, options is a structure with
% fields:
%   - **large** [true|{false}] if true, a computation explicitly using the
%   kronecker product is avoided.
%   - **multiply** [true|{false}]: if true, explicit omega matrices are
%   constructed and then multiplied to other matrices to sum the
%   permutations. Else, a functional form is used instead.
%
% Outputs
% --------
%
% - **res** [nd x nz^5]: output matrix
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if nargin<11
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

template=sparse(nd,nz^5);
res=template;

is_computable=@utils.cr.is_computable;

is_dv=is_computable(dv);
is_dvv=false;
is_dvvv=false;
is_dvvvv=false;
is_dvvvvv=false;

if is_dv
    is_dvv=is_computable(dvv);
    if is_dvv
        is_dvvv=is_computable(dvvv);
        if is_dvvv
            is_dvvvv=is_computable(dvvvv);
            if is_dvvvv
                is_dvvvvv=is_computable(dvvvvv);
            end
        end
    end
    if is_computable(vz)
        if is_computable(vzz)
            res=res+dvvv_vz_vzz_vzz();
            res=res+dvvvv_vz_vz_vz_vzz();
            if is_computable(vzzz)
                res=res+dvvv_vz_vz_vzzz();
                res=res+dvv_vzz_vzzz();
                if is_computable(vzzzz)
                    res=res+dvv_vz_vzzzz();
                    if is_computable(vzzzzz)
                        res=res+dv*vzzzzz;
                    end
                end
            end
        end
        res=res+dvvvvv_vz_vz_vz_vz_vz();
    end    
end

    function res=dvv_vzz_vzzz()
        res=template;
        if is_dvv
            if options.large
                res=utils.kronecker.A_times_kron_Q1_Qk(...
                    dvv,vzz,vzzz);
            else
                res=dvv*kron(vzz,vzzz);
            end
            
            if options.multiply
                res=res*cr.omega(nz,9);
            else
                res=utils.cr.dv_vz_omega(res,nz,9);
            end
        end
    end

    function res=dvv_vz_vzzzz()
        res=template;
        if is_dvv
            if options.large
                res=utils.kronecker.A_times_kron_Q1_Qk(...
                    dvv,vz,vzzzz);
            else
                res=dvv*kron(vz,vzzzz);
            end
            
            if options.multiply
                res=res*cr.omega(nz,8);
            else
                res=utils.cr.dv_vz_omega(res,nz,8);
            end
        end
    end

    function res=dvvv_vz_vzz_vzz()
        res=template;
        if is_dvvv
            if options.large
                res=utils.kronecker.A_times_kron_Q1_Qk(...
                    dvvv,vz,vzz,vzz);
            else
                res=dvvv*kronall(vz,vzz,vzz);
            end
            res=res*utils.cr.omega(nz,7);
        end
    end

    function res=dvvv_vz_vz_vzzz()
        res=template;
        if is_dvvv
            if options.large
                res=utils.kronecker.A_times_kron_Q1_Qk(...
                    dvvv,vz,vz,vzzz);
            else
                res=dvvv*kronall(vz,vz,vzzz);
            end
            if options.multiply
                res=res*cr.omega(nz,6);
            else
                res=utils.cr.dv_vz_omega(res,nz,6);
            end
        end
    end

    function res=dvvvv_vz_vz_vz_vzz()
        res=template;
        if is_dvvvv
            if options.large
                res=utils.kronecker.A_times_kron_Q1_Qk(...
                    dvvvv,vz,vz,vz,vzz);
            else
                res=dvvvv*kronall(vz,vz,vz,vzz);
            end
            if options.multiply
                res=res*cr.omega(nz,5);
            else
                res=utils.cr.dv_vz_omega(res,nz,5);
            end
        end
    end

    function res=dvvvvv_vz_vz_vz_vz_vz()
        res=template;
        if is_dvvvvv
            if options.large
                res=utils.kronecker.A_times_kron_Q1_Qk(...
                    dvvvvv,vz,vz,vz,vz,vz);
            else
                res=dvvvvv*utils.kronecker.kronall(vz,vz,vz,vz,vz);
            end
        end
    end
end