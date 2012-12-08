function [derivs,auxiliary,numberOfEquations,numberOfVariables,jac_toc]=...
    rise_derivatives(batch,validnames,wrt,Definitions)

% wrt can be entered in two ways:
% 1- a direct list of the variables to differentiate. e.g. wrt={'x','y','z'}
% 2- a two-column cell array with the first column containing the generic
% names of the variables and the second column the numeric indices. e.g.
% wrt={'x',1:10;'y',3:7}

if ischar(batch)
    batch=cellstr(batch);
end

if nargin>3
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
tic
[derivs,auxiliary]=sad_forward.jacobian(symbolic_batch,symb_list,with_respect_to);
derivs=trim_symbolic_equation(derivs);
jac_toc=toc();

derivs=analytical_symbolic_form(derivs,validnames,'analytic');
auxiliary=analytical_symbolic_form(auxiliary,validnames,'analytic');

numberOfEquations=numel(batch);
numberOfVariables=numel(with_respect_to);

