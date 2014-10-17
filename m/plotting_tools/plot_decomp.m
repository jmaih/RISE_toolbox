function h=plot_decomp(obj,nticks)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin<2
    nticks=[];
end

uall=double(obj);
upos=zeros(size(uall));
uneg=zeros(size(uall));
ipos=uall>=0;
ineg=uall<0;
upos(ipos)=uall(ipos);
uneg(ineg)=uall(ineg);
this_pos=ts(obj.start,upos,obj.varnames);
this_neg=ts(obj.start,uneg,obj.varnames);
this_all=ts(obj.start,sum(uall,2));
bar(this_pos,'nticks',nticks,'stack') %  area(this_pos)
hold on
bar(this_neg,'nticks',nticks,'stack') % area(this_neg)
hold on
plot(this_all,'k-','linewidth',2,'nticks',nticks)
axis tight;
hold off
if nargout>0
    h=gca;
end
