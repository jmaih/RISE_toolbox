function [x,x_id]=find_farthest(xvec,x0)
dd=abs(xvec-x0);
x_id=find(dd==max(dd),1,'first');
x=xvec(x_id);