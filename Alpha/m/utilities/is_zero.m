function flag=is_zero(a)
flag=(isnumeric(a)||islogical(a)) && all(a==0);
end
