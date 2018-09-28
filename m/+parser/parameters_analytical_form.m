function panal=parameters_analytical_form(pos)

psymb=parser.create_state_list('param',pos);

panal=parser.analytical_symbolic_form(psymb,'param','analytic');

end