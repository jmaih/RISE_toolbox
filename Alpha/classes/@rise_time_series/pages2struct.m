function this=pages2struct(this0)
vnames=this0.varnames;
if numel(unique(vnames))~=this0.NumberOfVariables
    error([mfilename,':: number of unique variable names different from number of columns of data matrix'])
end
this=struct();
datta=double(this0);
for ii=1:this0.NumberOfVariables
    this.(this0.varnames{ii})=rise_time_series(this0.start,permute(datta(:,ii,:),[1,3,2]));
end
end
