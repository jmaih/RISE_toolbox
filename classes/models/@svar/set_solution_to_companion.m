function [A,B,steady_state]=set_solution_to_companion(obj)
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
        error([mfilename,':: number of output arguments cannot exceed 1 when the object is empty'])
    end
    A=struct();
    return
end

if isempty(obj.solution)
    error('model has not been solved')
end

endo_nbr=obj.endogenous.number;
exo_nbr=sum(obj.exogenous.number);
reg_nbr=obj.markov_chains.regimes_number;

[A,~,B]=vartools.resolve(obj.solution,obj.nlags,reg_nbr);
cc=(obj.nlags-1)*endo_nbr;
nrows=obj.nlags*endo_nbr;
steady_state=cell(1,reg_nbr);
for ireg=1:reg_nbr
    steady_state{ireg}=repmat(obj.solution.ss{ireg},obj.nlags,1);
    % constant is at the end according to vartools.set_y_and_x
    %---------------------------------------------------------
    det_terms=A{ireg}(:,nrows+1:end-obj.constant);
    A{ireg}=[A{ireg}(:,1:nrows);
        eye(cc),zeros(cc,endo_nbr)];
    tmp=zeros(obj.endogenous.number,exo_nbr);
    tmp(:,obj.exogenous.is_observed)=det_terms;
    tmp(:,~obj.exogenous.is_observed)=B{ireg};
    B{ireg}=[tmp;
        zeros(cc,exo_nbr)];
end

end
