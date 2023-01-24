function plot_sens_in_data(grad,data)
%__________________________________________________________________________
% Script to plot ONLY the sensors in your data rather than all the channels
% in the grad structure.
%
% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)
%__________________________________________________________________________

grad_to_plot = [];
loc = ismember(grad.label,data.label);

grad_to_plot.label = grad.label(loc);
grad_to_plot.chanori = grad.chanori(loc,:);
grad_to_plot.chanpos = grad.chanpos(loc,:);
grad_to_plot.coilori = grad.coilori(loc,:);
grad_to_plot.coilpos = grad.coilpos(loc,:);

figure; ft_plot_sens(grad_to_plot,'facecolor','b',...
    'label','label','fontsize',5);
end
