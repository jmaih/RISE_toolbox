function obj=uplus(a)

if ~is_zero(a)
    obj=rise_sym(a);
else
    obj=rise_sym(0);
end
