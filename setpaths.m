function list=setpaths(flag)
if nargin<1
    flag=false;
end
filename=mfilename;
fullpath=which(filename);
loc=strfind(fullpath,filename);
tmp=genpath(fullpath(1:loc-2));

if ispc
    semicols=strfind(tmp,';');
elseif ismac
    semicols=strfind(tmp,':');
else
    error([mfilename,':: unknown system '])
end
previous=0;
collect={};
for ii=1:numel(semicols)
    current=semicols(ii);
    thepath=tmp(previous+1:current-1);
    if isempty(strfind(thepath,'.svn')) && ...
            isempty(strfind(lower(thepath),'test')) && ...
            isempty(strfind(lower(thepath),'junk'))
        if flag % then remove the paths
            rmpath(thepath)
        else
            collect=[collect,{thepath}];
            addpath(thepath)
        end
    end
    previous=current;
end
if nargout
	list=collect;
end
