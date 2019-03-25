%  INTERNAL FUNCTION: finds n tick positions in a vector of length "total"
% 
%  ::
% 
%    posticks=locate_ticks(total,n)
%    posticks=locate_ticks(total,n,start_at)
% 
%  Args:
% 
%     - **total** [integer]: length of the axis
%     - **n** [integer]: desired number of ticks
%     - **start_at** [integer|{1}]: location of first tick
% 
%  Returns:
%     :
% 
%     - **posticks** [vector]: location of the ticks
% 
%