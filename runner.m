%% EXAMPLE

close all; 
clear; clc;
%% Single file
sig = audioread("data/normal/85033_TV.wav");
labels = readtable("data/normal/85033_TV.tsv","FileType","delimitedtext","Delimiter","tab");

% props = project_run(sig);
props = get_features(sig, labels)';

%% Folder

%This folder only has 40058_TV and 85033_TV in this example
folder = 'data/normal/';
[props_n,labels_n] = run_for_folder(folder);

%This folder only has 61117_TV and 84802_PV in this example
folder = 'data/murmur/';
[props_m,labels_m] = run_for_folder(folder);

props = [props_n,props_m];
labels = [labels_n;labels_m];

%reading in ground truth HR-s (temp variable used to filter for these files)
hr_normal = readtable("data/HR_normal.csv");
hr_murmur = readtable("data/HR_murmur.csv");

temp = hr_normal.Signal=="40058_TV" | hr_normal.Signal=="85033_TV";
hr_n = table2array(hr_normal(temp,2));
temp = hr_murmur.Signal=="61117_TV" | hr_murmur.Signal=="84802_PV";
hr_m = table2array(hr_murmur(temp,2));
hrs = [hr_n;hr_m];

%setting ground truth array
pathology = [0,0,1,1];

%clearing variables
clear temp
clear props_true

%setting up the ground truth struct array
for k=1:length(hrs)
    temp.HR = hrs(k);
    temp.pathology = pathology(k);
    props_true(k) = temp;
end
%% Calculation
[hit_percent,miss_percent,multihit_percent,hrdiff_percent,ibsegdiff_percent,Se,Sp] = ...
    calc_score(props,props_true,labels);

avg_percent = mean([hit_percent;miss_percent;multihit_percent;hrdiff_percent;ibsegdiff_percent],2);
avg_hit_percent = avg_percent(1);
avg_miss_percent = avg_percent(2);
avg_multihit_percent = avg_percent(3);
avg_hrdiff_percent = avg_percent(4);
avg_ibsegdiff_percent = avg_percent(5);