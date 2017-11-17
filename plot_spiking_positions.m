%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s()
% Outputs:
% Inputs:
%   xN - normalized position x
%   yN - normalized position y
%   spikes_binned - spikes at corresponding xN and yN positions
%
% This function plots the spikes against the path taken as a function of
% normalized position.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_spiking_positions(xN, yN, XLocAtSpikes, YLocAtSpikes)
n = size(XLocAtSpikes,1);

for i = 1:n
    psp(i,xN,yN,XLocAtSpikes{i},YLocAtSpikes{i});
end

end

function psp(i,x, y, xs, ys)
figure('Name','Raw Data - Spikes vs Position');
plot(x,y,xs,ys,'r.');

% making it look pretty and readable
title([num2str(i) ': Spiking Activity as a Function of Position'])
xlabel('Normalized x position'); ylabel('Normalized y position');
legend('Path taken','Position at spike','location','best')
set(gca,'fontsize',14)

end