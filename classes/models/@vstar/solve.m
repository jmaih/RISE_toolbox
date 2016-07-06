function obj=solve(obj)
% SOLVE -- Builds the solution of a VSTAR object
%
% Syntax
% -------
% ::
%
%   obj=SOLVE(obj)
%
% Inputs
% -------
%
% - **obj** [vstar object]: model object
%
% Outputs
% --------
%
% - **obj** [vstar object]: model object with the state
% matrices in the field "solution"
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

obj.solution=get_parameters(obj);

end

function s=get_parameters(obj)

params=obj.parameter_values(obj.reordering_index);

s=struct();

p=obj.endogenous.number;

thresholds=obj.thresholds;

nthresh=numel(thresholds);

nregs=nthresh+1;

s.B=zeros(p,obj.num_regessors,nregs);

offset=0;

for ireg_=1:nregs
    
    s.B(:,:,ireg_)=reshape(pull_params(1:p*obj.num_regessors),p,obj.num_regessors);
    
end

nchol=1*p*(p+1)/2;

sigpars=pull_params(1:p);

corr_pars=pull_params(1:nchol-p);

s.C=eye(p);

for icol=1:p-1
    
    s.C(icol+1:end,icol)=corr_pars(1:p-icol);
    
    corr_pars(1:p-icol)=[];
    
end

s.C=s.C*diag(sigpars);

s.thresholds=cell(1,nthresh);

for ii=1:nthresh
    
    np=obj.thresholds(ii).np;
    
    myparams=pull_params(1:np);
    
    s.thresholds{ii}=struct('name',[obj.thresholds(ii).name,...
        '{',obj.thresholds(ii).lag,'}'],...
        'func',obj.thresholds(ii).func,...
        'g',myparams(1),'c',myparams(2:end));
    
end

    function pp=pull_params(span)
        
        stretch=offset+span;
        
        pp=params(stretch);
        
        offset=stretch(end);
        
    end

end
