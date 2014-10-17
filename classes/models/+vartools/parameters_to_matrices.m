function [links,parameter_values,all_param_names_vec]=parameters_to_matrices(param_template,...
    param_names,numberOfRegimes,debug)
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

if nargin<4
    debug=false;
    if nargin<3
        numberOfRegimes=1;
    end
end

% Make sure that the names are sorted
%------------------------------------
assert(isequal(param_names,sort(param_names)))

% link the parameters to the structural matrices : M(i,j)=p(k)
%-------------------------------------------------------------
ncell=size(param_template,2);
links=cell(1,ncell);

bigimage=param_template(2,:);
for ii=1:ncell
    themat=param_template{2,ii};
    header=param_template{1,ii};
    [jj,kk]=find(isnan(themat));
    if isempty(jj)
        continue
    end
    plist=cellstr(strcat(header,'_',int2str(jj),'_',int2str(kk)));
    ploc=locate_variables(plist,param_names);
    links{ii}=sub2ind(size(themat),jj,kk)+ploc*1i;
    image_themat=num2cell(themat);
    for jj=1:numel(links{ii})
        ll=links{ii}(jj);
        image_themat{real(ll)}=param_names{imag(ll)};
    end
    bigimage{ii}=image_themat;
    if debug
        disp(themat)
        disp(image_themat)
        keyboard
    end
end
bigimage=cellfun(@(x)replace_double(x),[bigimage{:}],'uniformOutput',false);
all_param_names_vec=bigimage(:);

nparams=numel(param_names);
parameter_values=nan(nparams,numberOfRegimes);

    function y=replace_double(x)
        y=x;
        if isa(y,'double')
            y='louisPergaudBahoyaMaih';
        end
    end
end