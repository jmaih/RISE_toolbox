function [obj,retcode]=solve(obj,varargin)%
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
    reg_nbr=obj.markov_chains.regimes_number;
    endo_nbr=obj.endogenous.number(end);
    steady_state=repmat({zeros(endo_nbr,1)},1,reg_nbr);
    for istate=1:reg_nbr
        for ilink=1:numel(links)
            mat=obj.param_template{2,ilink};
            mat(real(links{ilink}))=obj.parameter_values(imag(links{ilink}),istate);
                obj.solution.(obj.param_template{1,ilink}){istate}=mat;
        end
        if obj.constant
            % constant is at the end according to vartools.set_y_and_x 
            %---------------------------------------------------------
            constant=obj.solution.c{istate}(:,end);
            LEFT=eye(endo_nbr);
            if isfield(obj.solution,'a0')
                LEFT=obj.solution.a0{istate};
            end
            for ilag=1:obj.nlags
                lag_name=sprintf('a%0.0f',ilag);
                LEFT=LEFT-obj.solution.(lag_name){istate};
            end
            steady_state{istate}=LEFT\constant;
        end
    end
    obj.solution.ss=steady_state;
end

