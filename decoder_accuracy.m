% decoder_accuracy.m
%
% Script that uses decodes the position of a rat based using spiking activity
% from an ensemble of hippocampal place cells
%
% GROUP NAME: __________________________________
%
% Macauley Breault
% Created: 11-02-2017

clear; close all; clc;

% Load data
load('test.mat')

%% Decode position

% ~~~~~~~~~~~~~~~~~~~~ PLEASE ENTER YOUR DECODER FUNCTION HERE ~~~~~~~~~~~~~~~~~~~~
% 
% This function should have the following format
% 
% Inputs:
%     - spikes_binned = <Tx10> logical matrix of spiking activity in 1 ms time bins
%                         Assume it has the same format as train.mat
%                         but with a different number of time bins
% 
% Outpts:
%     - xN_decode = <Tx1> double vector of decoded x position of rat
%     - yN_decode = <Tx1> double vector of decoded y position of rat

[xN_decode, yN_decode] = your_decoder_function(spikes_binned);


%% Plot decoder results

figure
clf

% Plot
hold on
plot(xN,yN,'-b')
plot(xN_decode, yN_decode,'--r')
plot(cos(0:pi/50:2*pi),sin(0:pi/50:2*pi),'k-','LineWidth',2)
hold off

% Format
box('on')
axis('tight','square')

% Label
xlabel('X')
ylabel('Y')
legend('Actual','Decoded','Location','best')

%% Calculate accuracy using

% RMSE
accuracy = sqrt(sum(diag(cov(xN-xN_decode,yN-yN_decode)).^2));

disp(['Your decoder had an accuracy of ',num2str(accuracy)])
