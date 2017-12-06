%% load data
clr;

% glm_out-neuron_1.mat
% glm_out-neuron_2.mat
% 171130-glm_out-neuron_3.mat
% glm_out-neuron_4.mat
% 171130-glm_out-neuron_5.mat
% 171130-glm_out-neuron_6.mat
% 171130-glm_out-neuron_7.mat
% 171130-glm_out-neuron_8.mat
% 171130-glm_out-neuron_9.mat
% 171130-glm_out-neuron_10.mat
load('glm_out-neuron_1');

%% plot betas
figure(1); clf; hold on;
set(gcf,'units','points','position',[100,100,1000,600])

subplot(3,1,1);
errorbar(b3,2*stats3.se);
% xticks(1:length(b3));
xlim([0 length(b3)]);
xlabel('\beta number'); ylabel('\beta value');
set(gca,'FontSize',16)

% subplot(4,1,2);
% plot(b3);
% xlim([0 length(b3)]);
% xlabel('\beta number'); ylabel('\beta value');

%% plot p-values
subplot(3,1,2); hold on;
set(gca,'FontSize',16)
plot(stats3.p)
xlim([0 length(stats3.p)]);
% ylim([0 0.05])
ylabel('p-values'); xlabel('covariate number');

subplot(3,1,3); hold on;
plot(find(stats3.p<0.05),stats3.p(stats3.p<0.05), 'ro')
plot(find(stats3.p<0.05),stats3.p(stats3.p<0.05), 'b')
xlim([0 length(stats3.p)]);
ylabel('p-values'); xlabel('covariate number');
set(gca,'FontSize',16)

%%  find stretches of significant covariates
figure(2); clf; hold on;
set(gcf,'units','points','position',[100,100,1000,400])

subplot(1,2,1); hold on;
set(gca,'FontSize',16)
[L,n] = bwlabel(stats3.p < 0.05);
h = histogram(L,n);
xlabel('Stretch number'); ylabel('Length of stretch');
xlim([0 n]);

len = sort(h.Values);
longest = len(end);  % non-zero bins

subplot(1,2,2); hold on;
set(gca,'FontSize',16)
histogram(L,n);
xlabel('Stretch number'); ylabel('Length of stretch');
ylim([0 len(end-1)]); xlim([0 n]);

stretch_indexes = find(h.Values==longest)-1;
for i=1:numel(stretch_indexes)
    stretch_index = stretch_indexes(i);
    covar_index = find(L == stretch_index);

    for j=1:5
        covar_index = covar_index(covar_index ~= j);
    end

    hist_dep = covar_index - 5
end