function [Jac,i_index,myrefs]=Jacobian(func,wrt,ValidNames)
if nargin<3
    ValidNames={'x','y','param','def'};
    if nargin<2
        wrt=[];
    end
end

i_index=0;
myrefs={};
% myfunc='exp(x(1)+acos(x(2)*log(x(3)))+atan(x(1)*x(2)))';
if ischar(func)
    func=cellstr(func);
end

nrows=numel(func);
varlist=[];
for irow=1:nrows
    switch class(func{irow})
        case 'function_handle'
            func{irow}=func2str(func{irow});
        case 'char'
        otherwise
            error([mfilename,':: expecting a char or a function handle or a cell array of any of those'])
    end
    [func{irow},varlist_]=analytical_symbolic_form(func{irow},ValidNames,'symbolic');
    varlist=[varlist,varlist_];
end
varlist=unique(varlist);
if isempty(wrt)
    wrt=varlist;
end

[~,wrt]=analytical_symbolic_form(wrt,ValidNames,'symbolic');
ncols=numel(wrt);
%%
Jac=cell(nrows,ncols);
%%
func=cell2mat(strcat(varlist,','));
func=str2func(['@(',func(1:end-1),')',strout]);

for ii=1:numel(varlist)
    if ii==1
        varlist{ii}=sad_tree(varlist{ii},[],1);
    else
        varlist{ii}=sad_tree(varlist{ii});
    end
end

tree=func(varlist{:});
wrt=varlist;
myderivs=cell(2,numel(wrt));
for ii=1:numel(wrt)
    myderivs{1,ii}=wrt{ii}.name;
    myderivs{2,ii}=diff(tree,wrt{ii});
end
%% re_flag the tree after the derivatives have been computed
i_index=re_flag_tree(tree,i_index);
%%
myrefs=[myrefs;collect_references(tree)];
%% printing the derivatives
for ii=1:numel(wrt)
    myderivs{2,ii}=char(myderivs{2,ii});
myderivs{2,ii}=analytical_symbolic_form(myderivs{2,ii},ValidNames,'analytic');
end
% print(tree)