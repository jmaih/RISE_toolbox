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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

vnames=this0.varnames;

if numel(unique(vnames))~=this0.NumberOfVariables
    
    error([mfilename,':: number of unique variable names different ',...
        'from number of columns of data matrix'])
    
end

this=struct();

datta=this0.data;

for ii=1:this0.NumberOfVariables
    
    newdata=permute(datta(:,ii,:),[1,3,2]);
    
    if ii==1
        
        this.(this0.varnames{ii})=ts(this0.date_numbers,newdata,[],false,true);
        
    else
        
        this.(this0.varnames{ii})=reset_data(this.(this0.varnames{1}),newdata);
        
    end
    
end

end
