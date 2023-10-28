function [features,labels] = run_for_folder(folder, path)
%   Goes through each file and gets the feature of it
%   Requires that the wav and tsv files have the same name

    files = dir([folder '*.wav']);
    files = struct2table(files);

    for k=1:length(files.name)
        fname = files.name(k);
        fname = fname{:};
        sig = audioread([folder fname]);
        lab = readtable([folder fname(1:end-4) '.tsv'],"FileType","delimitedtext","Delimiter","tab");
        labels(k,:) = path;
        features(k,:) = get_features(sig);
    end

end