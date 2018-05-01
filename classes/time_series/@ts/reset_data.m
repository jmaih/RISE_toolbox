function this=reset_data(this,newdata,newnames)
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
%    - Keeps the start date and only changes the data. New data no longer
%    needs to be of the same size as the old one.
%
% Example:
%
%    See also:

nd=size(newdata,2);

comments=repmat({''},1,nd);

if nargin<3
    
    newnames=comments;
    
end

if ischar(newnames)
    
    newnames=cellstr(newnames);
    
end

n=numel(newnames);

if ~all(cellfun(@isempty,newnames,'uniformOutput',true)) && numel(unique(newnames))~=n
    
    error('some duplicated names')
    
end

if size(newdata,2)~=n
    
    error('number of names does not match number of new variables')
    
end

ts.check_size(newdata);

this.data=newdata;

this.varnames=newnames;

this.description=comments;

end
