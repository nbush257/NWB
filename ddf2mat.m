function ddf2mat(varargin)
%% function ddf2nwb([pname])
if length(varargin)==1
    pname = varargin{1};
    d_all = dir([pname '/*.ddf']);
else
    d_all = dir('*.ddf');
end

fnames = {d_all.name};
tag = {};
count = 1;
for ii = 1:length(fnames)
    
    idx = regexp(fnames(ii),'-\d{6}.ddf');
    iscalib = regexp(fnames(ii),'calib');
    if ~cellfun(@isempty,iscalib)
        continue
    else       
        tag{count} = fnames{ii}(1:idx{1}-1);
    end
    count = count + 1;
end
tag = unique(tag);

if length(varargin)==1
    p = pname;
else
    p = pwd;
end
%% make subdirectories for each cell/whisker
done_flag = false(size(tag));
for ii = 1:length(tag)
    target_dir = [p '/nwb_' tag{ii}];
    if isdir(target_dir);
        done_flag(ii)=true;
    else
        
    mkdir(target_dir);
    end
    
end
% Loop through each ddf
for ii = 1:length(tag)
    if done_flag(ii)
        warning('Experiment %s already converted, skipping',tag{ii});
        continue
    end
    
    fprintf('Converting:\t%s\n',tag{ii})
    % get a list of all ddfs from a particular cell/whisker
    recordings = dir([p '\' tag{ii} '*.ddf']);
    
    % concatenate all the raw data to vizualize only

    for jj = 1:length(recordings)
        outfname = [target_dir '\' recordings(jj).name(1:end-4) '_tempConvert.mat'];
        data = importDDF_vg([p '\' recordings(jj).name]);
        save(outfname,'data');
       
    end

end

