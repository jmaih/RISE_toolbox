%  PROGRESS_MONITOR :  
% 
% 
% 
% 
%  OUTPUT
%  -------
%  -**D** [parallel.pool.DataQueue] : parallel.pool.DataQueue object to be
%  used directly or to construct a handle to monitor progress
%  * direct usage : send(D,[ichain,ismpl])
%  * handle usage : control=@(ismpl)send(D,[ichain,ismpl]);
% 
% 
%  MAIN (CORE) INPUTS
%  ------------------
%  -**distribution** [vector] : number of iterations to be run by each
%    worker
% 
%  VARARGIN INPUTS
%  ----------------
%  -**nSegments** [integer|{10}] : total number of segments
% 
%  -**ProgressSymbols** [cell array|{{'-', '\\', '|', '/'}}] : flipping
%  symbols during progress
% 
%  -**CompletedSymbol** [char|{'~'}] : symbols indicating completion
% 
%  -**TodoSymbol** [char|{'-'}] : symbols indicating incomplete elements
% 
%  -**DelimiterSymbols** [cell array|{{'[', ']'}}] : symbols delineating the
%  print out on screen
% 
%  Examples
%  --------
% 
%  D=progress_monitor(nLoops*ones(1,nchain));
% 
%  D=progress_monitor(nLoops*ones(1,nchain),'DelimiterSymbols', {'(', ')'}, ...
%      'ProgressSymbols', {'o', '0', 'O', '0'}, ...
%      'CompletedSymbol', '>', 'TodoSymbol', '<');
% 
%  D=progress_monitor(nLoops*ones(1,nchain),'DelimiterSymbols', {'[', ']'}, ...
%      'ProgressSymbols', {'o', '0', 'O', '0'}, ...
%      'CompletedSymbol', 'd', 'TodoSymbol', 'b');
%