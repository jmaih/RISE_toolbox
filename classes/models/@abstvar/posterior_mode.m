function dd=posterior_mode(self)

dd=[self.estim_.links.estimList(:),num2cell(self.estim_.estim_param)];

dd=cell2struct(dd(:,2),dd(:,1),1);


end