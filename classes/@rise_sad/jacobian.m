function [Jac,... % jacobian in string form with auxiliary terms
    dd,... % jacobian in rise_sad form. Readily usable for further differentiation
    references,... % definition of the auxiliary terms
    tank ... % containers.Map object for the different operations in the objectives
    ]=jacobian(objectives,... set of objective functions (rise_sad, string or function handles)
    varnames,... variable names entering the objectives
    wrt,... list of variables to differentiate with respect to
    references,...
    tank ...
    )

if nargin<5
    tank=[];
    if nargin<4
        references=[];
        if nargin<3
            wrt=[];
            if nargin<2
                error('At least two arguments should be provided')
            end
        end
    end
end
vlist='';
args=cell2rise_sad(varnames);
if isempty(wrt)
    wrt=args;
end
%% create tree
switch class(objectives)
    case 'rise_sad'
        tree=objectives(:);
    case {'char','cell'}
        if ischar(objectives)
            objectives=cellstr(objectives);
        end
        tree=rise_sad.empty(0);
        n=0;
        for ic=1:numel(objectives)
            n=n+1;
            switch class(objectives{ic})
                case 'rise_sad'
                    tree(n,1)=objectives{ic};
                case 'char'
                    if isempty(vlist)
                        vlist=rise_sad2cellstr(varnames);
                        vlist=unique(vlist);
                        vlist=cell2mat(strcat(vlist(:)',','));
                        vlist=vlist(1:end-1);
                    end
                    objectives{ic}=str2func(['@(',vlist,')',objectives{ic}]);
                    tree(n,1)=objectives{ic}(args{:});
                case 'function_handle'
                    tree(n,1)=objectives{ic}(args{:});
                otherwise
                    error(['Unsupported class ',class(objectives{ic})])
            end
        end
    case 'function_handle'
        tree=objectives(args{:});
    otherwise
        error(['Unsupported class ',class(objectives)])
end
%% put wrt in the right form
wrt=cell2rise_sad(wrt);
%% compute derivatives... and register the number of calls to each node in the tree
[dd,references,tank]=diff(tree,wrt,references,tank);
%% now the derivatives can be printed
Jac=cell(size(dd));
for irow=1:size(dd,1)
    for jcol=1:size(dd,2)
        Jac{irow,jcol}=char(dd(irow,jcol));
    end
end

end

function dd=cell2rise_sad(item)
if ischar(item)
    item=cellstr(item);
end
dd=cell(size(item));
for it=1:numel(item)
    dd{it}=rise_sad(item{it});
end
end

function item=rise_sad2cellstr(item)
switch class(item)
    case 'rise_sad'
        item={item.name};
    case 'cell'
        for it=1:numel(item)
            if ~ischar(item{it})
                item{it}=rise_sad2cellstr(item{it});
            end
        end
    case 'char'
        item=cellstr(item);
    otherwise
        error(['Unsupported class ',class(item)])
end
end
