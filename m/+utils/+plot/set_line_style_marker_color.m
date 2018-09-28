function set_line_style_marker_color(h,clrs,lnwdt)

% h is a vector of handles to the function plot
% clrs is a matrix rgb
% lnwdt is a numeric for line width

if nargin<3
    
    lnwdt=[];
    
    if nargin<2
        
        clrs=[];
        
    end
    
end

lineStyle={'-','--',':','-.'};

marker={'o','+','*','.','x','s','d','^','v','>','<','p','h'};

if isempty(clrs)
    
    clrs=gray(15); % OK: gray, hsv,bone,copper
    
end

if isempty(lnwdt),lnwdt=1.5; end

cliter=0;

mliter=0;

lsliter=0;

nl=numel(lineStyle);

nm=numel(marker);

nc=size(clrs,1);

for ii=1:numel(h)
    
    cliter=cliter+1; mliter=mliter+1; lsliter=lsliter+1;
    
    mk=marker{mliter};
    
    ls=lineStyle{lsliter};
    
    cl=clrs(cliter,:);
    
    set(h(ii),'lineStyle',ls,'marker',mk,'color',cl,'linewidth',lnwdt,...
        'MarkerSize',1.5);
    
    if mliter==nm,mliter=0;end
    
    if cliter==nc,cliter=0;end
    
    if lsliter==nl,lsliter=0;end
    
end

end