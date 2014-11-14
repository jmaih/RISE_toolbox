function retcode=check_derivatives(obj,varargin)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also:


if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    retcode=struct();
    return
end
DerivativesTypes={'symbolic','numerical','automatic'};
nder=numel(DerivativesTypes);

nobj=numel(obj);
retcode_=nan(1,nobj);
for imod=1:nobj
    retcode_(imod)=recheck_derivatives();
end
if nargout
    retcode=retcode_;
end


    function rcode=recheck_derivatives()
        obj(imod)=set(obj(imod),varargin{:});
        solve_order=obj(imod).options.solve_order;
        
        evaluated_objects=cell(1,nder);
        structural_matrices=cell(1,nder);
        rcode=0;
        ider=0;
        while ~rcode && ider<nder
            ider=ider+1;
            tic
            [evaluated_objects{ider},rcode,structural_matrices{ider}]=...
                solve(obj(imod),'solve_derivatives_type',DerivativesTypes{ider});
            tt=toc;
            if rcode
                disp([DerivativesTypes{ider},...
                    ' derivatives or solution failed for model(',int2str(imod),...
                    ') with error message '])
                decipher(rcode);
            else
                disp(['model(',int2str(imod),') computation of derivatives and solution under ',...
                    DerivativesTypes{ider},' derivatives::',num2str(tt),' seconds'])
            end
        end
        
        if ~rcode
            dv='d';
            Tz='T';
            main=['model(',int2str(imod),')::'];
            for io=1:solve_order
                dv=[dv,'v']; %#ok<*AGROW>
                Tz=[Tz,'z'];
                order_=int2str(io);
                for ider=1:nder-1
                    if io>1 && strcmp(DerivativesTypes{ider},'numerical')
                        fprintf('numerical derivatives not available for order %0.0f\n',io)
                    else
                        disp([main,DerivativesTypes{ider},' derivatives(order=',order_,'): discrep of derivatives relative to automatic'])
                        test=do_test(dv,1);
                        disp([dv,'::',test])
                    end
                    
                    if obj(imod).is_optimal_policy_model && io==1
                        disp([main,DerivativesTypes{ider},' derivatives(order=',order_,'): discrep of policy weights relative to automatic'])
                        test=do_test('weights',3);
                        disp(['Weights ::',test])
                    end
                    
                    if io>1 && strcmp(DerivativesTypes{ider},'numerical')
                        fprintf('solution using numerical derivatives not available for order %0.0f\n',io)
                    else
                        disp([main,DerivativesTypes{ider},' derivatives(order=',order_,'): discrep or Solution relative to automatic'])
                        test=do_test(Tz,2);
                        disp([Tz,'::',test])
                    end
                end
            end
        end
        function discr=do_test(vx,flag)
            switch flag
                case 1
                    discr=cell2mat(structural_matrices{ider}.(vx))-...
                        cell2mat(structural_matrices{nder}.(vx));
                case 3
                    discr=cell2mat(structural_matrices{ider}.planner.(vx))-...
                        cell2mat(structural_matrices{nder}.planner.(vx));
                case 2
                    discr=cell2mat(evaluated_objects{ider}.solution.(vx))-...
                        cell2mat(evaluated_objects{nder}.solution.(vx));
            end
            discr=num2str(max(abs(discr(:))));
        end
    end

end

