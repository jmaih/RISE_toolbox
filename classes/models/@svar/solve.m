function [obj,retcode]=solve(obj,varargin)%

if isempty(obj)
    obj=struct();
    return
end

if ~isempty(varargin)
    obj=set(obj,varargin{:});
end
nobj=numel(obj);
if nobj>1
    retcode=nan(1,nobj);
    for iobj=1:nobj
        [obj(iobj),retcode(iobj)]=solve(obj(iobj));
    end
    return
end

[obj.solution.transition_matrices,retcode]= ...
    utils.code.evaluate_transition_matrices(obj.routines.transition_matrix,obj.parameter_values(:,1));
if ~retcode
    links=obj.param_to_mat_links;
    for istate=1:obj.markov_chains.regimes_number
        for ilink=1:numel(links)
            mat=obj.param_template{2,ilink};
            mat(real(links{ilink}))=obj.parameter_values(imag(links{ilink}),istate);
                obj.solution.(obj.param_template{1,ilink}){istate}=mat;
        end
    end 
    if isa(obj,'rfvar')
        obj=structural_form(obj);
    end
end

