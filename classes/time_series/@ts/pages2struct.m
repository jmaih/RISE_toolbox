function this=pages2struct(this0)
% Turns multivariable ts object into a struct with time series object
%
% ::
%
%    output = pages2struct(input)
%
% Args:
%    input (ts object): ts object to turn into struct form of data
%
% Returns:
%    :
%
%    - **output** (struct): a struct with
%
%       - fieldname: variable names of input
%       - value: ts object corresponding to the variable (with data and description)
%

vnames=this0.varnames;

if numel(unique(vnames))~=this0.NumberOfVariables
    
    error([mfilename,':: number of unique variable names different ',...
        'from number of columns of data matrix'])
    
end

this=struct();

datta=this0.data;

for ii=1:this0.NumberOfVariables
    
    newdata=permute(datta(:,ii,:),[1,3,2]);
    
    % take the description in cell form
    description=this0.description(ii);
    
    if ii==1
        
        this.(this0.varnames{ii})=ts(this0.date_numbers,newdata,[],...
            description,true);
        
    else
        
        this.(this0.varnames{ii})=set(this.(this0.varnames{1}),'data',newdata);
                
    end
    
end

end
