function O=occurrence_map(batch,type,nv,collapse)
% When O is 3-dimensional, 
% - the first page are the lags,
% - the second page are the current,
% - the third page are the leads

if nargin<4
    
    collapse=false;
    
end

n=numel(batch);

C=regexp(batch,['\<',type,'\((?<vname>\d+)\)(?<lag>(\{(\+|-)?\d+\})?)'],'names');

O=false(n,nv,3);

for ii=1:n
    
    Ci=C{ii};
    
    for jj=1:numel(Ci)
        
        Cij=Ci(jj);
        
        vpos=str2double(Cij.vname);
        
        lag=2;
        
        if ~isempty(Cij.lag)
            
            lag=lag+str2double(Cij.lag(2:end-1));
            
        end
        
        if vpos>nv
            
            warning(['expansion occurs from ',int2str(nv),...
                ' to ',int2str(max(Cii))])
            
        end
        
        % never mind the repetitions in Cii %     Cii=unique(Cii);
        O(ii,vpos,lag)=true;
        
    end
    
end

if collapse
    
    O=logical(sum(O,3));
    
end

end

%{
function [O]=occurrence_map(batch,type,nv)

n=numel(batch);

% remove time
C=regexprep(batch,'\{(\+|-)?\d+\}','');

% extract type(n)
C=regexp(C,['\<',type,'\(\d+\)(\{(\+|-)?\d+\})?'],'match');

reconvert=@(x)cellfun(@(z)str2double(z),x,'uniformOutput',true);

O=false(n,nv,3);

for ii=1:n
    
    Cii=strrep(C{ii},[type,'('],'');
    
    Cii=reconvert(strrep(Cii,')',''));
    
    if ~isempty(Cii)
        
        if max(Cii)>nv
            warning(['expansion occurs from ',int2str(nv),...
                ' to ',int2str(max(Cii))])
        end
        
        % never mind the repetitions in Cii %     Cii=unique(Cii);
        O(ii,Cii)=true;
        
    end
    
end

end
%}