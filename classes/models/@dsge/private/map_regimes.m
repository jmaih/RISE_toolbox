%--- help for map_regimes ---
%
%  INTERNAL FUNCTION: map regimes from small ones to big ones and vice versa.
%  Useful when dealing with loose commitment.
% 
%  ::
% 
%    map_=map_regimes(obj)
%    map_=map_regimes(obj,big2small)
% 
%  Args:
%     obj (dsge | rise): model object
%     big2small ({true} | false): decides the direction of the mapping
% 
%  Returns:
%     :
% 
%     - **map_** [1xh|1xbigh]: vector assigning the regimes
% 
%
%    Other functions named map_regimes
%
%       dsge/map_regimes
%