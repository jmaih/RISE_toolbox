function this=pages2struct(this0)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

vnames=this0.varnames;
if numel(unique(vnames))~=this0.NumberOfVariables
    error([mfilename,':: number of unique variable names different from number of columns of data matrix'])
end
this=struct();
datta=this0.data;
for ii=1:this0.NumberOfVariables
    this.(this0.varnames{ii})=ts(this0.date_numbers,permute(datta(:,ii,:),[1,3,2]));
end
end
