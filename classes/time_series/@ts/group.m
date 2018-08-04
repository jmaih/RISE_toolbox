function g=group(this,varargin)
% GROUP -- groups contributions
%
% ::
%
%
%   g=GROUP(this,{group1,{v11,v12,...}},...,{groupn,{vn1,vn2,...}})
%
% Args:
%
%    - **this** [ts]: time series of multiple variables
%
%    - **varargin** []: sequences of cell arrays with the following possible
%    formats:
%      - {'supply','Ep'}
%      - {'demand',{'Ey','Er'}}
%
% Returns:
%    :
%
%    - **g** [ts]: new time series with grouped contributions
%
% Note:
%
% Example:
%
%    See also:


n=length(varargin);
discard=false(1,n);
for ii=1:length(varargin)
    if ischar(varargin{ii})
        discard(ii)=true;
        is_pop=strcmp(varargin{ii},'pop');
        if ~is_pop
            error('To discard the rest, the option is "pop"')
        end
    end
end
varargin=varargin(~discard);
n=length(varargin);

data=double(this);

start=this.start;

oldnames=this.varnames;

is_taken=false(1,this.NumberOfVariables);

nobs=this.NumberOfObservations;

newnames=cell(1,n);

newdata=nan(nobs,n);

for ii=1:n
    
    if numel(varargin{ii})~=2
        
        error('each group must be a two-element cell: e.g. {''supply'',{''E1'',''E2'',''E3''}}')
        
    end
    
    [newnames{ii},pos]=set_one_group();
    
    newdata(:,ii)=sum(data(:,pos),2);
    
    is_taken(pos)=true;
end

% add the unclassified
%---------------------
newdata=[newdata,data(:,~is_taken)];

newnames=[newnames,oldnames(:,~is_taken)];

g=ts(start,newdata,newnames);


    function [header,pos]=set_one_group()
        
        header=varargin{ii}{1};
        
        items=varargin{ii}{2};
        
        if ischar(items)
            
            items=cellstr(items);
            
        end
        
        if ~iscellstr(items)
            
            error('group elements must be char or cellstr')
            
        end
        
        if ~ischar(header)
            
            error('the first element of each cell must be a char')
            
        end
        
        pos=locate_variables(items,oldnames);
        
        bad=is_taken(pos);
        
        badnames=oldnames(bad);
        
        if ~isempty(badnames)
            
            disp(badnames)
            
            error('the variables above belong to several groups')
            
        end
        
    end

end