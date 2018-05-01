function params=isnan(obj)
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

if isempty(obj)
    
    params=cell(0,4);
    
    return
    
end

nobj=numel(obj);

if nobj>1
    
    params=cell(1,nobj);
    
    for iobj=1:numel(obj)
        
        params{iobj}=isnan(obj(iobj));
        
        disp(['====== model ',int2str(iobj),' ======'])
        
        disp(params{iobj})
        
    end
    
    return
    
end

pp=obj.parameter_values;

regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));

regimes_number=obj.markov_chains.regimes_number;

chain_names=obj.markov_chains.chain_names;

params={};

for ipar=1:sum(obj.parameters.number)
    
    nanlocs=isnan(pp(ipar,:));
    
    if any(nanlocs)
        
        chain_loc=obj.parameters.governing_chain(ipar);
        
        gov_chain=chain_names{chain_loc};
        
        pname=obj.parameters.name{ipar};
        
        if strcmp(gov_chain,'const')
            
            params=[params,pname]; %#ok<*AGROW>
        
        else
            % put nanlocs and states side by side and trim to have a
            % clearer picutre.
            states_nans=[regimes(:,chain_loc),nanlocs(:)];
            
            discard=false(regimes_number,1);
            
            states=[];
            
            for irow=1:regimes_number
                
                if ismember(states_nans(irow,1),states)
                    
                    discard(irow)=true;
                
                else
                    
                    states=[states,states_nans(irow,1)];
                
                end
                
            end
            
            states_nans=states_nans(~discard,:);
            
            for irow=1:size(states_nans,1)
                
                if states_nans(irow,2)
                    
                    params=[params,[pname,'(',gov_chain,',',sprintf('%0.0f',states_nans(irow,1)),')']];
                
                end
                
            end
            
        end
        
    end
    
end
