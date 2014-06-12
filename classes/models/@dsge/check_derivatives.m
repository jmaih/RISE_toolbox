function retcode=check_derivatives(obj,varargin)

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
retcode=nan(1,nobj);
for imod=1:nobj
    retcode(imod)=recheck_derivatives();
end


    function rcode=recheck_derivatives()
        obj(imod)=set(obj(imod),varargin{:});
        
        evaluated_objects=cell(1,nder);
        rcode=0;
        ider=0;
        while ~rcode && ider<nder
            ider=ider+1;
            tic
            [evaluated_objects{ider},rcode]=...
                evaluate(obj(imod),'solve_derivatives_type',DerivativesTypes{ider});
            tt=toc;
            if rcode
                disp([DerivativesTypes{ider},' derivatives failed for model(',int2str(imod),') with error message '])
                utils.error.decipher(rcode);
            else
                disp(['model(',int2str(imod),') evaluation under ',DerivativesTypes{ider},' derivatives::',num2str(tt),' seconds'])
            end
        end
        
        if ~rcode
            VariablesToCheck={'Aplus','A0','Aminus','B'};
            if obj(imod).is_optimal_policy_model
                VariablesToCheck=[VariablesToCheck,'W'];
            end
            [evaluated_objects{nder},rcode]=solve(evaluated_objects{nder});
            FurtherVariablesToCheck={'T','R'};
            for ider=1:nder-1
                disp(['model(',int2str(imod),')::',DerivativesTypes{ider},' derivatives relative to automatic derivatives'])
                for ivar=1:numel(VariablesToCheck)
                    test=evaluated_objects{ider}.(VariablesToCheck{ivar})-...
                        evaluated_objects{nder}.(VariablesToCheck{ivar});
                    disp([VariablesToCheck{ivar},'::',num2str(max(max(max(abs(test)))))])
                end
                if ~rcode
                    evaluated_objects{ider}=solve(evaluated_objects{ider});
                    for ivar=1:numel(FurtherVariablesToCheck)
                        test=evaluated_objects{ider}.(FurtherVariablesToCheck{ivar})-...
                            evaluated_objects{nder}.(FurtherVariablesToCheck{ivar});
                        disp([FurtherVariablesToCheck{ivar},'::',num2str(max(max(max(abs(test)))))])
                    end
                end
            end
        end
        
    end
end

