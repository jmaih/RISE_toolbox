function flag=is_one(a)
flag=(isnumeric(a)||islogical(a)) && all(a==1);
end
