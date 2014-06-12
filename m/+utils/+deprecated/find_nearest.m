function n=find_nearest(a,b)
dist=utils.optim.distance(a,b);
n=find(dist==min(dist),1,'first');