function c=char(x,long)
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

if nargin<2
    long=false;
end

if isnumeric(x)
    n=numel(x.func);
    c=num2str(x.func);
    if n>1
        c=['[',c,']'];
    end
elseif isempty(x.args)
    c=x.func;
else
    args=x.args;
    for iarg=1:numel(args)
        args{iarg}=char(args{iarg},long);
    end
    the_func=x.func;
    % use vectorized operation directly
    %----------------------------------
    if any(strcmp(the_func,{'mpower','mrdivide','mtimes'}))% <-- ismember(the_func,{'mpower','mrdivide','mtimes'})
        the_func=the_func(2:end);
    end
    if strcmp(the_func,'plus')
        c=[args{1},'+',args{2}];
    elseif long||any(strcmp(the_func,{'abs','acos','acosh','asin','asinh','atan','atanh','cos',...
            'cosh','cot','erf','exp','log','log10','sign','sin','sinh',...
            'tan','if_elseif','if_then_else','normcdf','normpdf','max','min','sqrt'}))
        c=cell2mat(strcat(args,','));
        c=[the_func,'(',c(1:end-1),')'];
% % % % % % 		if any(strcmp(the_func,{'if_elseif','if_then_else'}))
% % % % % % 	        c=['utils.functional_programming.',c];
% % % % % % 		end
    else
        switch the_func
            case 'uplus'
                c=args{1};
            case 'uminus'
                c=['-(',args{1},')'];
            case 'and'
                c=['(',args{1},' & ',args{2},')'];
            case 'eq'
                c=['(',args{1},')==(',args{2},')'];
            case 'ge'
                c=['(',args{1},')>=(',args{2},')'];
            case 'gt'
                c=['(',args{1},')>(',args{2},')'];
            case 'le'
                c=['(',args{1},')<=(',args{2},')'];
            case 'lt'
                c=['(',args{1},')<(',args{2},')'];
            case 'minus'
                c=[args{1},'-(',args{2},')'];
            case 'power'
                % use the vectorized form directly
                c=['(',args{1},').^(',args{2},')'];
            case 'rdivide'
                % use the vectorized form directly
                c=['(',args{1},')./(',args{2},')'];
            case 'times'
                % use the vectorized form directly
                c=['(',args{1},').*(',args{2},')'];
            case 'ne'
                c=['(',args{1},')~=(',args{2},')'];
            case 'or'
                c=['(',args{1},')|(',args{2},')'];
            otherwise
                error(['unknown function :: ',the_func])
        end
    end
end