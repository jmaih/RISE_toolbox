function check_shock_id(id,nx)
id=id(:);
if any(id>nx)||...
        any(id<=0)||...
        any(floor(id)~=id)
    msg=['id should be an integer between 1 and ',int2str(nx)];
    error(msg)
end
end
