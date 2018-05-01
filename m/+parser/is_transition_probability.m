function [istp,diagonal,chain_name,max_state]=is_transition_probability(name)
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

if iscellstr(name)
    
    n=numel(name);
    
    if n==1
        
        name=name{1};
        
    else
        
        istp=false(1,n);
        
        diagonal=false(1,n);
        
        chain_name=cell(1,n);
        
        max_state=zeros(1,n);
        
        for ii=1:n
            
            [istp(ii),diagonal(ii),chain_name{ii},max_state(ii)]=...
                parser.is_transition_probability(name{ii});
            
        end
        
        return
        
    end
    
end

istp=false;

diagonal=false;

chain_name='';

max_state=0;

if isletter(name(1))

    underscore=strfind(name,'_');
    
    if ~isempty(underscore)&& numel(underscore)==3
    
        if size(name,2)>=8 && strcmp('tp',name(underscore(1)+1:underscore(2)-1))
        
            first=name(underscore(2)+1:underscore(3)-1);
            
            second=name(underscore(3)+1:end);
            
            if ~all(isletter(first)) && ~all(isletter(second))
            
                first=str2double(first);
                
                second=str2double(second);
                
                if ~isnan(first) && ~isnan(second)
                
                    istp=true;
                    
                    chain_name=name(1:underscore(1)-1);
                    
                    max_state=max(first,second);
                    
                    if second==first
                    
                        diagonal=true;
                    
                    end
                    
                end
                
            end
            
        end
        
    end
    
end