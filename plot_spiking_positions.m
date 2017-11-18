%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_spiking_positions(xN,yN,XLocAtSpikes,YLocAtSpikes)
% Outputs: none, displays n figures
% Inputs:
%   xN - normalized positions x (of any # of neurons)
%   yN - normalized positions y (of any # of neurons)
%   XLocAtSpikes - positions xN at spikes (of any # of neurons)
%   YLocAtSpikes - positions yN at spikes (of any # of neurons)
%
% This function plots the spikes and the path taken as a function of
% normalized position for all neurons given. It loops through each neuron
% and calls the function psp(i,xN,yN,xspikes_i,yspikes_i) to plot
% corresponding figures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_spiking_positions(xN,yN,XLocAtSpikes,YLocAtSpikes)
n = size(XLocAtSpikes,1);

% loop through each neuron
for i = 1:n
    psp(i,xN,yN,XLocAtSpikes{i},YLocAtSpikes{i});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psp(i,xN,yN,xspikes_i,yspikes_i)
% Outputs: none, displays figure
% Inputs:
%   i - neuron ID #
%   xN - normalized positions x
%   yN - normalized positions y
%   xspikes_i - positions xN at spikes (XLocAtSpikes{i})
%   yspikes_i - positions yN at spikes (YLocAtSpikes{i})
%
% This function plots the spikes against the path taken as a function of
% normalized position. The figure is labeled according to neuron ID #. This
% function is called by
% plot_spiking_positions(xN,yN,XLocAtSpikes,YLocAtSpikes) for n number of
% neurons in the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psp(i,xN,yN,xspikes_i,yspikes_i)
figure('Name','Raw Data - Spikes and Path vs Position');
plot(xN,yN,xspikes_i,yspikes_i,'r.');

% making it look pretty and readable
title([num2str(i) ': Spiking Activity as a Function of Position'])
xlabel('Normalized x position'); ylabel('Normalized y position');
% legend('Path taken','Position at spike','location','best')
set(gca,'fontsize',14)

end