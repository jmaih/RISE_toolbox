function retcode=check_derivatives(obj,varargin)
% check_derivatives - compares the derivatives and the solutions from various differentiation techniques
%
% ::
%
%
%   check_derivatives(obj)
%   retcode=check_derivatives(obj)
%
% Args:
%
%    - **obj** [rise|dsge]: model object or vectors of model objects
%
% Returns:
%    :
%
%    - **retcode** [numeric]: 0 if no problem is encountered during the
%      comparisons. Else the meaning of recode can be found by running
%      decipher(retcode)
%
% Note:
%
%    - The derivatives computed are 'automatic', 'symbolic' or 'numerical'
%
%    - The comparisons are done relative to automatic derivatives, which are
%      assumed to be the most accurate.
%
% Example:
%
%    See also:


if isempty(obj)
    
    if nargout>1
        
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
        
    end
    
    retcode=cell(0,4);
    
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
        
        tt=nan(1,nder);
            
        main=['model(',int2str(imod),')::'];
        
        fprintf('\n ******************* %s%s ******************* \n',main,obj(imod).filename)
        
        while ~rcode && ider<nder
            
            ider=ider+1;
            
            tic
            
            [evaluated_objects{ider},rcode,structural_matrices{ider}]=...
                solve(obj(imod),'solve_derivatives_type',DerivativesTypes{ider});
            
            tt(ider)=toc;
            
            if rcode
                
                disp([DerivativesTypes{ider},...
                    ' derivatives or solution failed for model(',int2str(imod),...
                    ') with error message '])
                decipher(rcode);
                
            else
                
                disp(['model(',int2str(imod),') derivatives and solution under ',...
                    DerivativesTypes{ider},' derivatives::',num2str(tt(ider)),' seconds'])
                
            end
            
        end
        
        if ~rcode
            
            dv='d';
            
            Tz='T';
            
            for io=1:solve_order
                
                dv=[dv,'v']; %#ok<*AGROW>
                
                Tz=[Tz,'z'];
                
                order_=int2str(io);
                
                for ider=1:nder-1
                    
                    if io>1 && strcmp(DerivativesTypes{ider},'numerical')
                        
                        fprintf('numerical derivatives and solution not available for order %0.0f\n',io)
                        
                    else
                        
                        fprintf('\n%s %s derivatives(order = %s)\n',main,DerivativesTypes{ider},order_)
                        
                        test=do_test(dv,1);
                        
                        fprintf('* discrep of derivatives relative to automatic : %s :: %s \n',dv,test)
                        
                        if (obj(imod).is_optimal_policy_model||...
                                obj(imod).is_optimal_simple_rule_model) && io==1
                            
                            test=do_test('weights',3);
                            
                            fprintf('* discrep of policy weights relative to automatic : Weights :: %s \n',test)
                            
                        end
                        
                        test=do_test(Tz,2);
                        
                        fprintf('* discrep or Solution relative to automatic : %s :: %s \n',Tz,test)
                        
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

