function this=reset_data(this,newdata,newnames)
% Reset the data
%
% ::
%
%    db = reset_data(db,newdata);
%    db = reset_data(db,newdata,newnames);
%
% Args:
%    - db (ts object): original time series object
%    - newdata (matrix): new data to replace the data of the original time series
%    - newnames: (optional) new variable names to replace with
%
% Returns:
%    :
%    db (ts object): time series object with the data (and variable names) replaced
%
% Note:
%    - Keeps the start date and only changes the data. New data does not needs to be of the same size as the old one.
%

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
