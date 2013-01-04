function retcode=check_derivatives(obj,varargin)

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    retcode=struct();
    return
end

obj=set_options(obj,varargin{:});

DerivativesTypes={'symbolic','numerical','automatic'};
nder=numel(DerivativesTypes);
evaluated_objects=cell(1,nder);
retcode=0;
ider=0;
while ~retcode && ider<nder
    ider=ider+1;
    tic
    [evaluated_objects{ider},retcode]=...
        evaluate(obj,'derivatives',DerivativesTypes{ider});
    tt=toc;
    if retcode
        disp([DerivativesTypes{ider},' derivatives failed with error message '])
        decipher_error(retcode);
    else
        disp(['model evaluation under ',DerivativesTypes{ider},' derivatives::',num2str(tt),' seconds'])
    end
end

if ~retcode
    VariablesToCheck={'Aplus','A0','Aminus','B'};
    if obj.is_optimal_policy_model
        VariablesToCheck=[VariablesToCheck,'W'];
    end
    [evaluated_objects{nder},retcode]=solve(evaluated_objects{nder});
    FurtherVariablesToCheck={'T','R'};
    for ider=1:nder-1
        disp([DerivativesTypes{ider},' derivatives relative to automatic derivatives'])
        for ivar=1:numel(VariablesToCheck)
            test=evaluated_objects{ider}.(VariablesToCheck{ivar})-...
                evaluated_objects{nder}.(VariablesToCheck{ivar});
            disp([VariablesToCheck{ivar},'::',num2str(max(max(max(abs(test)))))])
        end
        if ~retcode
            evaluated_objects{ider}=solve(evaluated_objects{ider});
            for ivar=1:numel(FurtherVariablesToCheck)
                test=evaluated_objects{ider}.(FurtherVariablesToCheck{ivar})-...
                    evaluated_objects{nder}.(FurtherVariablesToCheck{ivar});
                disp([FurtherVariablesToCheck{ivar},'::',num2str(max(max(max(abs(test)))))])
            end
        end
    end
end


