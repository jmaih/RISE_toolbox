function g=group(this,groups)
% Groups contributions
%
% ::
%
%   g=group(this,{group1,{v11,v12,...}},...,{groupn,{vn1,vn2,...}})
%
% Args:
%
%    this (ts): time series of multiple variables
%
% - groups : [structure|cell array |{empty}] grouping of shocks in the decomposition.
%   By default, the shocks are not grouped. The syntax is of the form
%   {group1,{v11,v12,...},...,groupn,{vn1,vn2,...}}. The shocks that are
%   not listed are put in a special group called "others". The "others"
%   group does not include the effect of initial conditions.
%   e.g. p=struct();
%        p.demand={'Ey','Er'};
%        p.supply={'Ep'};
%   e.g. p={'demand',{'Ey','Er'},'supply',{'Ep'}};
%
% Returns:
%    :
%
%    - **g** [ts]: new time series with grouped contributions
%

oldnames=this.varnames;

[groups,pos_groups]=utils.miscellaneous.check_groups(oldnames,groups);

data=double(this);

start=this.start;

nobs=this.NumberOfObservations;

newnames=fieldnames(groups);

n=numel(newnames);

newdata=nan(nobs,n);

for ii=1:n
    
    pos=pos_groups{ii};
    
    newdata(:,ii)=sum(data(:,pos),2);
    
end

g=ts(start,newdata,newnames);

end