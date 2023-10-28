close all; 
clear; clc;

%% running classification on whole folder 
folder = 'data/normal/';
files = dir([folder '*.wav']);
files = struct2table(files);
for k=1:length(files.name)
    fname = files.name(k);
    fname = fname{:};
    sig = audioread([folder fname]);
    features(k,:) = get_features(sig);
    props(k) = get_props(sig, features(k,:), 0);
end

%% running classification on  single signal
fname = "data/normal/85033_TV.wav";
sig = audioread(fname);
features = get_features(sig)';
props = get_props(sig, features, 1);


%%
close all; 
clear; clc;
%reading in ground truth HR-s (temp variable used to filter for these files)
folder = 'data/normal/';
files = dir([folder '*.wav']);
files = struct2table(files);

for k=1:length(files.name)
        fname = files.name(k);
        fname = fname{:};
        sig = audioread([folder fname]);
        lab = readtable([folder fname(1:end-4) '.tsv'],"FileType","delimitedtext","Delimiter","tab");
        lab = table2cell(lab);
        labels_n{k,:} = lab;
        features(k,:) = get_features(sig);
        props_n(k) = get_props(sig, features(k,:), 0);
        fprintf('%s\n',fname);
        pathology(1,k) = 0;

 end

folder = 'data/murmur/';
files = dir([folder '*.wav']);
files = struct2table(files);

for k=1:length(files.name)
        fname = files.name(k);
        fname = fname{:};
        sig = audioread([folder fname]);
        lab = readtable([folder fname(1:end-4) '.tsv'],"FileType","delimitedtext","Delimiter","tab");
        lab = table2cell(lab);
        labels_m{k,:} = lab;
        features(k,:) = get_features(sig);
        props_m(k) = get_props(sig, features(k,:), 0);
        fprintf('%s\n',fname);
        pathology(1,50+k) = 0;
 end

props = [props_n,props_m];
labels = [labels_n;labels_m];

hr_normal = readtable("data/HR_normal.csv");
hr_murmur = readtable("data/HR_murmur.csv");

hr_n = table2array(hr_normal(:,2));
hr_m = table2array(hr_murmur(:,2));
hrs = [hr_n;hr_m];

%setting ground truth array

%clearing variables
clear temp
clear props_true

%setting up the ground truth struct array
for k=1:length(hrs)
    temp.HR = hrs(k);
    temp.pathology = pathology(k);
    props_true(k) = temp;
end
[hit_percent,miss_percent,multihit_percent,hrdiff_percent,ibsegdiff_percent,Se,Sp] = ...
    calc_score(props,props_true,labels);

avg_percent = mean([hit_percent;miss_percent;multihit_percent;hrdiff_percent;ibsegdiff_percent],2);
avg_hit_percent = avg_percent(1);
avg_miss_percent = avg_percent(2);
avg_multihit_percent = avg_percent(3);
avg_hrdiff_percent = avg_percent(4);
avg_ibsegdiff_percent = avg_percent(5);
