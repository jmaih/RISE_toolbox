function d=transpose(db)
d=main_frame(db,false);

d=permute(d,[2,1,3]);
end