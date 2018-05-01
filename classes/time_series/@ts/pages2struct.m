function this=pages2struct(this0)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

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
        
        this.(this0.varnames{ii})=reset_data(this.(this0.varnames{1}),newdata);
        
        this.(this0.varnames{ii}).description=description;
        
    end
    
end

end
