function []=dim_relevant(freq,nber_datesfunc,figs)
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

if nargin < 3
    
    figs=[];
    
    if nargin < 2
        
        nber_datesfunc = [];
                
        if nargin<1
            
            freq=[];
            
        end
        
    end
    
end

if isempty(figs)
    
    figs=gcf;
    
end

nfigs=numel(figs);

if nfigs>1
    
    for ifig=1:nfigs
        
        dim_relevant(freq,nber_datesfunc,figs(ifig))
        
    end
    
    return
    
end

if isempty(nber_datesfunc)
    
    [start,finish]= nber_dates(freq);
    
else
        
    nber_datesfunc=fcnchk(nber_datesfunc);
    
    [start,finish]= nber_datesfunc(freq); %#ok<RHSFN>
    
end

if ischar(start)
    
    start= cellstr(start);
    
end

if ischar(finish)
    
    finish = cellstr(finish);
    
end
% shadenber.m
%
%  Routine to shade the nberdates in a figure

set_bounds()

colorstr=[159 182 205]/256;

dim(start,finish,colorstr,figs);

    function set_bounds()
        
        figure(figs)
        
        curax = axis;
        
        start_dn = date2serial(start);
        
        finish_dn = date2serial(finish);
        
        if start_dn(1)<curax(1)
            
            begin=find(start_dn>=curax(1),1,'first');  % First recession to include;
            
        else
            
            begin=1;
            
        end
        
        if finish_dn(end)>curax(2)
            
            last=find(finish_dn<=curax(2),1,'last');  % last recession to include;
            
        else
            
            last=numel(finish_dn);
            
        end
        
        start=start_dn(begin:last);
        
        finish=finish_dn(begin:last);
        
    end

end

