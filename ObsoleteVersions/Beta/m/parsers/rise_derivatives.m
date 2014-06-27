function [finalOutput,numberOfEquations,numberOfVariables,jac_toc]=...
    rise_derivatives(batch,validnames,wrt,Definitions,order)

% wrt can be entered in two ways:
% 1- a direct list of the variables to differentiate. e.g. wrt={'x','y','z'}
% 2- a two-column cell array with the first column containing the generic
% names of the variables and the second column the numeric indices. e.g.
% wrt={'x',1:10;'y',3:7}

if nargin<5
    order=[];
    if nargin<4
        Definitions=[];
    end
end

if isempty(order)
    order=1;
end
if order>2
    error('Forbidden to compute derivatives beyond 2')
end
if ischar(batch)
    batch=cellstr(batch);
end

if ~isempty(Definitions)
    batch=replace_definitions(batch,Definitions);
end

% symbolic form
symbolic_batch=analytical_symbolic_form(batch,validnames,'symbolic');

% list of symbols
symb_list=collect_symbolic_list(symbolic_batch,strcat(validnames,'_'));

if ischar(wrt)
    wrt=cellstr(wrt);
end

reprocess=size(wrt,2)==2 && isnumeric(wrt{1,2});
if reprocess
    with_respect_to={};
    for ii=1:size(wrt,1)
        digits=wrt{ii,2};
        xx=wrt{ii,1};
        for id=1:numel(digits)
            with_respect_to=[with_respect_to,[xx,'_',int2str(digits(id))]]; %#ok<AGROW>
        end
    end
else
    with_respect_to=wrt;
end

symb_list=union(symb_list,with_respect_to);

numberOfEquations=numel(symbolic_batch);
numberOfVariables=numel(with_respect_to);

tic
if order==1
    test=false;
    % vectorized form
    if test
        vectform=false;
    else
        vectform=true;
    end
    [~,~,~,finalOutput]=sad_reverse.jacobian(symbolic_batch,symb_list,with_respect_to,vectform);
    argouts={'Jac_'};
elseif order==2
    % non-vectorized form
    [~,~,~,finalOutput]=sad_reverse.hessian(symbolic_batch(1),symb_list,with_respect_to);
    argouts={'Hess_','Jac_'};
end
jac_toc=toc();

finalOutput=analytical_symbolic_form(finalOutput,validnames,'analytic');

finalOutput=struct('code',cell2mat(finalOutput(:)'),'argins',...
    {validnames},'argouts',{argouts});


