function a=alpha_probability(f_theta,f0)
r=exp(f_theta-f0);
a=min(1,r);