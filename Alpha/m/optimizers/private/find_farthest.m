function n=find_farthest(a,b)
dist=distance(a,b);
n=find(dist==max(dist),1,'first');