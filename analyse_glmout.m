%% load data
clr;

load('glm_out-neuron_3.mat');

%% plot betas
figure(1); clf; hold on;
set(gcf,'units','points','position',[100,100,1000,400])
subplot(1,2,1);
errorbar(b3,2*stats3.se);
xticks(1:length(b3));
xlim([0 length(b3)+1]);
xlabel('\beta number'); ylabel('\beta value');
subplot(1,2,2);
plot(b3);
xlabel('\beta number'); ylabel('\beta value');

%% plot p-values
figure(2); clf; hold on;
set(gcf,'units','points','position',[100,100,1000,400])
subplot(1,2,1); hold on;
plot(stats3.p)
% ylim([0 0.05])
ylabel('p-values'); xlabel('covariate number');
subplot(1,2,2); hold on;
plot(stats3.p(stats3.p<0.05), 'r.')
ylabel('p-values'); xlabel('covariate number');

%%  find stretches of significant covariates
figure(3); clf; hold on;
[L,n] = bwlabel(stats3.p < 0.05);
h = histogram(L,n);
xlabel('Stretch number'); ylabel('Length of stretch');

len = sort(h.Values);
longest = len(end-1);

ylim([0 longest])

stretch_indexes = find(h.Values==longest)-1;
for i=1:numel(stretch_indexes)
    stretch_index = stretch_indexes(i);
    covar_index = find(L == stretch_index);

    for i=1:5
        covar_index = covar_index(covar_index ~= i);
    end

    hist_dep = covar_index - 5
end