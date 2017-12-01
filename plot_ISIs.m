% plot_ISIs.m
% This function makes a histogram of interspike intervals (ISIs) for each
% Neuron in the dataset. It takes a matrix of spike trains (spikes) and the
% max ISI to include in the histogram (threshold) as inputs and returns a 
% cell array containing all the ISIs for each neuron.

function ISIs = plot_ISIs(spikes,threshold,binwidth)

n = size(spikes,2);

ISIs = cell(n,1);
figure('units','normalized','outerposition',[0 0.035 1 0.92]); 

for i = 1:n,
    
    ISIs{i} = diff(find(spikes(:,i)));
    pISI = ISIs{i}(ISIs{i}<=threshold); %Throw out ISIs that are too large (for plotting only)
    subplot(2,5,i); grid on; hold on;
    histogram(pISI,'BinWidth',binwidth);
    xlabel('Interspike Interval (ms)');
    ylabel('Number');
    title(['Neuron ' num2str(i)]);
    set(gca,'fontsize',20)

end
set(gca,'fontsize',20)
% suptitle('Distribution of Interspike Intervals');
end


