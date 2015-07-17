function index(self)
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

if ~isempty(self.date_numbers)
    fprintf(1,'start_date: %s end_date: %s frequency: %s\n',self.start,self.finish,self.frequency);
end
fprintf(1,'# variables: %0.0f # observations: %0.0f # pages: %0.0f \n',self.NumberOfVariables,self.NumberOfObservations,self.NumberOfPages);

end