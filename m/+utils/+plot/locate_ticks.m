function posticks=locate_ticks(total,n,start_at)
% locate_ticks -- finds n tick positions in a vector of length "total"
%
% ::
%
%
%   posticks=locate_ticks(total,n)
%
%   posticks=locate_ticks(total,n,start_at)
%
% Args:
%
%    - **total** [integer]: length of the axis
%
%    - **n** [integer]: desired number of ticks
%
%    - **start_at** [integer|{1}]: location of first tick
%
% Returns:
%    :
%
%    - **posticks** [vector]: location of the ticks
%
% Note:
%
% Example:
%
%    See also:

if nargin<3
    
    start_at=1;
    
end

incmnt=0;

while 1
    
    incmnt=incmnt+1;
    
    posticks=start_at:incmnt:total;
    
    if numel(posticks)<=n
        
        break
        
    end
    
end

end