% plot_ISIs.m
% This function makes a histogram of interspike intervals (ISIs) for each
% Neuron in the dataset. It takes a matrix of spike trains (spikes) and the
% max ISI to include in the histogram (threshold) as inputs and returns a 
% cell array containing all the ISIs for each neuron.

function ISIs = plot_ISIs(spikes,threshold)
n = size(spikes,2);

ISIs = cell(n,1);
figure;
for i = 1:n,
    ISIs{i} = diff(find(spikes(:,i)));
    pISI = ISIs{i}(ISIs{i}<=threshold); %Throw out ISIs that are too large (for plotting only)
    subplot(2,5,i); grid on; hold on;
    histogram(pISI,'BinWidth',10);
    xlabel('Interspike Interval (ms)');
    ylabel('Number');
    title(strcat('Distribution of Interspike Intervals for neuron: ',num2str(i)));    
end

end

