function batch=read_identification_restrictions(self,restrictions,shocknames)

endonames=self.endogenous;

[n,p]=size(restrictions);

if ~iscell(restrictions)|| p~=2
    
    error('restrictions should be a cell array with two columns')
    
end

batch=struct('vbl',{},'shk',{},'lag',{},'type',{},'rhs',{});

for ii=1:n
    
    lhs=restrictions{ii,1};
    
    rhs=restrictions{ii,2};
    
    [batch(ii).vbl,batch(ii).shk,batch(ii).lag,...
        batch(ii).type,batch(ii).rhs]=read_restriction(lhs,rhs);

end

    function [vbl,shk,lag,type,rhs]=read_restriction(lhs,rhs)
        
        at=find(lhs=='@');
        
        lag=0;
        
        if isempty(at)
            
            error(['error in restriction ',lhs])
            
        end
        
        vbl=lhs(1:at-1);
        
        shk=lhs(at+1:end);
        
        lp=find(vbl=='{');
        
        rp=find(vbl=='}');
        
        if ~((isempty(lp) && isempty(rp))||(~isempty(lp) && ~isempty(rp)))
            
            error('opening or closing parenthesis not matched')
            
        end
        
        if ~isempty(lp)
            
            % there could be many lags a:b, [a,b,c,d]
            
            lag=eval(vbl(lp+1:rp-1)); %<-- lag=str2double(vbl(lp+1:rp-1));
            
            if any(isnan(lag))
                
                error(['wrong lag specification in restriction "',lhs,'"'])
                
            end
            
            vbl=vbl(1:lp-1);
            
        end
        
        vbl=locate_variables(vbl,endonames);
        
        is_ls=all(isstrprop(shk,'digit'));
        
        if is_ls
            
            type='ls';
            
            shk=str2double(shk);
            
        else
            
            shk=locate_variables(shk,shocknames);
            
        end
        
        if isnumeric(rhs)
            
            if isscalar(rhs) && rhs==0
                
                if ~is_ls
                    
                    type='z';
                    
                end
                
            elseif numel(rhs)==2
                
                rhs=sort(rhs);
                
                type='m';
                
            else
                
                error(['unable to identifiy rhs of restriction "',lhs,'"'])
                
            end
            
        elseif ischar(rhs)
            
            comparison=strcmp(rhs,{'+','-'});
            
            if any(comparison)
                
                type='s';
                
                rhs=if_then_else(comparison(1),1,-1);
                
            else
                
                error(['unable to identifiy rhs of restriction "',lhs,'"'])
                
            end
            
        else
            
            error(['unable to identifiy rhs of restriction "',lhs,'"'])
            
        end
        
    end

end