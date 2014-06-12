function param_template=set_structural_matrices_structure(model_class,n,nlags,nx,iblx) % block exogenous indices
% the nan terms are the ones that are to be estimated
% the diagonal terms of the a0 matrix are moved to SIG so as to have a way
% of normalizing the switching of volatility. A natural question is whether
% a0 and SIG should switch simultaneously. Not switching them
% simultaneously implies a lower computational burden.

iblx=relogicalize(iblx,'exogenous block');
endo=~iblx;

% istat=relogicalize(istat,'stationary');
is_svar=strcmp(model_class,'svar');
% [a0,a1,...,ap,OMEGA,SIG,det]
%-----------------------------
param_template=cell(3,nlags+is_svar+2+nx>0);

offset=0;
% a0
%---
if is_svar
    param_template{1,1}='a0';
    param_template{2,1}=diag(ones(n,1));
    param_template{2,1}(param_template{2,1}==0)=nan;
    param_template{2,1}(iblx,endo)=0;
    param_template{3,1}=[n,n];
    offset=offset+1;
end

% a1, a2,...,ap
%--------------
for ilag=1:nlags
    param_template{1,offset+1}=sprintf('a%0.0f',ilag);
    atmp=nan(n);
    atmp(iblx,endo)=0;
    param_template{2,offset+1}=atmp;
    param_template{3,offset+1}=[n,n];
    offset=offset+1;
end

% OMEGA
%------
param_template{1,offset+1}='omg';
param_template{2,offset+1}=eye(n);
if ~is_svar
    param_template{2,offset+1}=param_template{2,offset+1}+tril(nan(n),-1);
end
param_template{3,offset+1}=[n,n];
offset=offset+1;

% SIG
%----
param_template{1,offset+1}='sig';
param_template{2,offset+1}=diag(nan(n,1));
param_template{3,offset+1}=[n,n];
offset=offset+1;

% deterministic terms
%--------------------
if nx
    param_template{1,offset+1}='c';
    param_template{2,offset+1}=nan(n,nx);
    param_template{3,offset+1}=[n,nx];
end

    function logic=relogicalize(logic,type)
        if isempty(logic)
            logic=false(1,n);
        end
        if ~islogical(logic)
            if any(logic>n)||any(logic<1)||any(fix(logic)~=logic)
                error(['wrong specification of ',type,' indices'])
            end
            tmp=false(1,n);
            tmp(logic)=true;
            logic=tmp;
        end
    end
end