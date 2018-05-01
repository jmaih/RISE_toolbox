function IDs=locate_variables(Variables,State,silent,NewAlgo)
% locate_variables locate variables in an array
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

if nargin<4
    
    NewAlgo=true;
    
    if nargin<3
        
        silent=false;
        
    end
    
end

if ischar(State)
    
    State=cellstr(State);
    
end

if isempty(Variables)
    
    nvar=0;
    
elseif ischar(Variables)
    
    nvar=size(Variables,1);
    
    Variables=cellstr(Variables);
    
elseif iscellstr(Variables)
    
    nvar=numel(Variables);
    
end

IDs=nan(nvar,1);

for j=1:nvar
    
    vj=Variables{j};
    % I remove spaces in the variables just to make sure... I hope this
    % will not come back to haunt me. (March 09, 2012)
%     vj(isspace(vj))=[]; 
%     % this came back to haunt me on May 17, 2012 and so I removed it.
    if NewAlgo
        
        vj_id=find(strcmp(vj,State));
        
    else
        
        vj_id=strmatch(vj,State,'exact'); 
        
    end
    
    if isempty(vj_id)
        
        if silent
            
            continue
            
        else
            
            error([mfilename,':: variable ',vj,' not found'])
            
        end
        
    elseif numel(vj_id)>1
        
        if silent
            
            continue
            
        else
            
            error([mfilename,':: variable ',vj,' declared more than once'])
        
        end
        
    end
    
    IDs(j)=vj_id;
    
end

% Replacements:
% Slow: STRMATCH(String, CellString)
% Fast: FIND(STRNCMP(String, CellString, length(String)))
%
% Slow: STRMATCH(String, CellString, 'exact')
% Fast: And FIND(STRCMP(String, CellString))

% tic
% for i=1:1e+5
%     IDs1=locate_variables(M_.endo_names,M_.endo_names(oo_.dr.order_var,:),true,true);
% end
% toc
% Elapsed time is 12.084246 seconds.
% tic
% for i=1:1e+5
%     IDs2=locate_variables(M_.endo_names,M_.endo_names(oo_.dr.order_var,:),true,false);
% end
% toc
% Elapsed time is 48.853927 seconds.


