function [pop,order]=sort_population(pop)
violstrength=[pop.violstrength];
fs=[pop.f];
order=[];
uv=unique(violstrength);
for iv=1:numel(uv)
    locs=find(violstrength==uv(iv));
    [~,tags]=sort(fs(locs));
    order=[order,locs(tags)]; %#ok<AGROW>
end
pop=pop(order);
end
