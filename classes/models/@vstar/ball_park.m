function B=ball_park(obj)

variables_locations_in_data=obj.variables_locations_in_data;

endo_data=obj.data(variables_locations_in_data.endo_id,:);

exo_data=obj.data(variables_locations_in_data.det_id,:);

[y,X]=vartools.set_y_and_x(endo_data,exo_data,obj.nlags,obj.constant);
        
B=y/X; % <- y*X'*inv(X*X')

end