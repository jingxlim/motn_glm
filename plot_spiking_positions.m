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
function plot_spiking_positions(xN,yN,XLocAtSpikes,YLocAtSpikes,varargin)

p = inputParser;

defaultType = 'subplot';
validType = {'subplot','separate'};
checkType = @(x) any(validatestring(x,validType));

addRequired(p,'xN');
addRequired(p,'yN');
addRequired(p,'XLocAtSpikes');
addRequired(p,'YLocAtSpikes');
addOptional(p,'type',defaultType,checkType);

parse(p,xN,yN,XLocAtSpikes,YLocAtSpikes,varargin{:})

if strcmp(p.Results.type,'separate')

    for i=1:numel(XLocAtSpikes)
        psp(i,xN,yN,XLocAtSpikes{i},YLocAtSpikes{i});
    end

elseif strcmp(p.Results.type,'subplot')

    figure('Name','Raw Data - Spikes and Path vs Position',...
        'units','normalized','outerposition',[0 0.035 1 0.92]); hold on;
    for i=1:numel(XLocAtSpikes)
        
        subplot(numel(XLocAtSpikes)/5, 5, i); hold on;
        plot(xN,yN,XLocAtSpikes{i},YLocAtSpikes{i},'r.');
        
        % make pretty
        title(['Neuron ' num2str(i)])
        xlabel('Normalized x position'); ylabel('Normalized y position');
        pbaspect([1 1 1]);
    end
    suptitle('Spiking activity as a function of position')
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
title(['Neuron: ' num2str(i) ': Spiking Activity as a Function of Position'])
xlabel('Normalized x position'); ylabel('Normalized y position');
% legend('Path taken','Position at spike','location','best')
set(gca,'fontsize',14)

end