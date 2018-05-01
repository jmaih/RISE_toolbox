function out=index(self)
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

sp=cell(0,1);

if ~isempty(self.date_numbers)
    
    sp=[sp,...
        sprintf('start_date: %s end_date: %s frequency: %s',...
        self.start,self.finish,self.frequency)];
    
end

sp=[sp,...
    sprintf('# variables: %0.0f # observations: %0.0f # pages: %0.0f',...
    self.NumberOfVariables,self.NumberOfObservations,self.NumberOfPages)];

if nargout
    
    out=sp;
    
else
    
    for ii=1:numel(sp)
        
        fprintf(1,'%s \n',sp{ii});
        
    end
    
end

end