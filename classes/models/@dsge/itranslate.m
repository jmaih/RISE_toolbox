function outList=itranslate(obj,inList,order)
% translate -- Translates comprehensible atoms (model variables) into RISE codes
%
% ::
%
%   outList=translate(obj,inList)
%   outList=translate(obj,inList,order)
%
% Args:
%
%    obj (rise | dsge): scalar model object.
%
%    inList (char | cellstr): List of atoms to translate e.g. C,
%       X{+1}, lambda_x, EPS_A,
%
%    order (integer|0|{1}): If order>1 the returned list of atoms is a list
%       of kroneckers in which the number of elements in each items is the
%       order of the kronecker. This is useful for instance to understand
%       what combination of variables makes up a column in a
%       differentiation. If order=0, the translation is with respect to the
%       static, rather than the dynamic model
%
% Returns:
%    :
%
%    - **outList** [cellstr]: List of translated atoms
%
% Example:
%    :
%
%    list=get(obj,'endo_list')
%    outList=itranslate(obj,list)
%    outList=itranslate(obj,list,0) % contemporaneous or steady state
%
%    list=get(obj,'exo_list')
%    outList=itranslate(obj,list)
%
%    list=get(obj,'param_list')
%    outList=itranslate(obj,list)
%
%    list=get(obj,'def_list')
%    outList=itranslate(obj,list)
%
% See also : dsge/translate
%    :
%
%

if isempty(obj)
    
    outList=cell(0,4);
    
    return
    
end

if nargin<3
    
    order=1;
    
end

if ischar(inList)
    
    inList=cellstr(inList);
    
end

w=yxpd();

outList=inList;

incidence=obj.lead_lag_incidence.after_solve;

expr='\<[a-zA-Z]+\w*\>';

expr=['(?<atom>',expr,')(\(|\{)?(?<digit>(\+|\-)?\d+)?(\)|\})?'];

p=regexp(inList,expr,'names');

p=[p{:}];

n=numel(p);

for ii=1:n
    
    t=p(ii).atom;
    
    try
        
        descr=w.(t);
        
    catch
        
        error([t,' not recognized as an atom in the model'])
        
    end
    
    loc=descr{1};
    
    teipe=descr{2};
    
    if strcmp(teipe,'y')
        
        v=str2double(p(ii).digit);
        
        if order==0
            
            if ~(isnan(v)||v==0)
                
                error('no leads or lags at order 0')
                
            end
            
        else
            
            if isnan(v)||v==0 % contemporaneous by default
                
                loc=incidence(loc,2);
                
            elseif v==1
                
                loc=incidence(loc,1);
                
            elseif v==-1
                
                loc=incidence(loc,3);
                
            end
            
        end
        
    end
    
    outList{ii}=sprintf('%s(%d)',teipe,loc);
    
end


    function [w,xList,pList,yList,dList]=yxpd()
        
        w=struct();
        
        imap=struct('xList','x',...
            'pList','param',...
            'yList','y',...
            'dList','def');
        
        xList=get(obj,'exo_list');
        
        pList=get(obj,'par_list');
        
        yList=get(obj,'endo_list');
        
        dList=get(obj,'def_list');
        
        absorb(xList,pList,yList,dList)
        
        function absorb(varargin)
            
            for jjj=1:length(varargin)
                
                list=varargin{jjj};
                
                ipn=inputname(jjj);
                
                for iii=1:numel(list)
                    
                    v=list{iii};
                    
                    w.(v)={iii,imap.(ipn)};
                    
                end
                
            end
            
        end
        
    end

end