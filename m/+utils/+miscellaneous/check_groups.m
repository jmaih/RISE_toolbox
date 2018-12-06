function [groups,positions]=check_groups(sList,groups)
% check_groups : checks the adequacy of the partitions for any
% decomposition
%
% Syntax
% -------
% ::
%
%   [groups,positions]=check_groups(sList)
%   [groups,positions]=check_groups(sList,groups)
% Inputs
% -------
%
% - sList : [char|cell array] list of the variables(contributing factors);
%
% - groups : [structure|cell array |{empty}] grouping of shocks in the decomposition.
%   By default, the shocks are not grouped. The syntax is of the form
%   {group1,{v11,v12,...},...,groupn,{vn1,vn2,...}}. The shocks that are
%   not listed are put in a special group called "others". The "others"
%   group does not include the effect of initial conditions.
%   e.g. p=struct();
%        p.demand={'Ey','Er'};
%        p.supply={'Ep'};
%   e.g. p={'demand',{'Ey','Er'},'supply',{'Ep'}};
%
% Outputs
% --------
%
% - groups : [struct] structure with the groupings. If the original groups
%   is empty, the decomposition is element by element. If there are
%   variables not included in the decomposition, a field called "others"
%   that gathers all of the unclassified variables.
%
% - positions : [cell array] each cell contains a vector indexes/locations
%   of the partition 
%
% Examples
% ---------
%
% See also: dsge/historical_decomposition

if nargin<2
    
    groups=[];
    
end

if ischar(sList)
    
    sList=cellstr(sList);
    
end

% horizontalize in case groups is empty
sList=sList(:).';

if isempty(groups)
    % groups of one variable
    %-------------------------    
    groups=cell2struct(sList,sList,2);
    
end

nshocks=numel(sList);

isLocated=false(1,nshocks);

structuralize()

group_names=fieldnames(groups);

ngroups=numel(group_names);

positions=cell(1,ngroups);

check_adequacy()

others=sList(~isLocated);

if ~isempty(others)
    
    groups.others=others;
    
    positions=[positions,{find(~isLocated)}];
    
end

    function check_adequacy()
               
        for ig=1:ngroups
            
            if strcmpi(group_names{ig},'others')
                
                error('No group in the decomposition can be called others')
                
            end
            
            vList=groups.(group_names{ig});
            
            if ischar(vList)
                
                vList={vList};
                
            end
            
            pos=locate_variables(vList,sList,true);
            
            bad=isnan(pos);
            
            if any(bad)
                
                disp(vList(bad))
                
                error('The names above are not in the main list')
                
            end
            
            already_assigned=isLocated(pos);
            
            if any(already_assigned)
                
                disp(vList(already_assigned))
                
                error('The names above belong to at least 2 groups')
                
            end
            
            isLocated(pos)=true;
            
            positions{ig}=pos;
            
        end
        
    end

    function structuralize()
        
        if isstruct(groups)
            
            return
            
        end
        
        two_n=numel(groups);
        
        if rem(two_n,2)
            
            error('The number of elements in the group should be even')
            
        end
        
        groups=reshape(groups,[],2);
        
        newgroup=struct();
        
        for ii=1:two_n/2
            
            gname=groups{ii,1};
            
            if ~isvarname(gname)
                
                error('"%s" is not a valid variable name',parser.any2str(gname))
                
            end
            
            vnames=groups{ii,2};
            
            newgroup.(gname)=vnames;
            
        end
        
    end

end