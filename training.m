% This script creates the feature vectors to be used for training

close all; 
clear; clc;

% getting the feature vectors of normal signals
folder = 'data/normal/';
[features_n,labels_n] = run_for_folder(folder,0);

% getting the feature vectors of murmur signals
folder = 'data/murmur/';
[features_m,labels_m] = run_for_folder(folder,1);

% adding the feature vectors together
features = [features_n;features_m];
labels = [labels_n;labels_m];

% adding label 
data = array2table([zeros(100,1),features]);
data.Properties.VariableNames(1) = {'Label'};
data.Label = labels;
