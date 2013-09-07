function n=find_nearest(a,b)
dist=distance(a,b);
n=find(dist==min(dist),1,'first');