function db=head(db,n)
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

if nargin<2
    n=min(5,db.NumberOfObservations);
end

db=ts(db.date_numbers(1:n),db.data(1:n,:,:),db.varnames);

end