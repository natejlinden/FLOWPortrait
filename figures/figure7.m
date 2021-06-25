% ORIGINAL CODE FROM Linden et al. 2020 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)

% This script recreates all the panes in Figure 7
% The script does not compute FTLE on every segment of data, it load pre-computed FTLE fields. The data sizes and computational times are too large to include here.

% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

clear all; close all; clc

% Plots FLOW portraits accross all subjects and developmental days
% Computes and plots quantification metrics 

addpath('../flow_portraits/'); 


dir_list={
    '051118_P1_Ani2_Ca_00002',
    '120817_animal1_P1_CaData00001',
    '050218_P2_Ani2_Ca_00001',
    '050218_P2_Ani2_Ca_00002',
    '042318_P3_Ani1_00002',
    '042318_P3_Ani1_00004',
    '042318_P3_Ani1_00005',
    '042318_P3_Ani1_00007',
    '042518_P5_Ani1_00001',
    '042518_P5_Ani3_Ca_00001_67_98281',
    '120617_animal4_P7_CaData00019',
    '120617_animal4_P7_CaData00020',
    '121017_Animal_3_P7_CaData00001',
    '121017_Animal_3_P7_CaData00002',
    '121017_Animal_3_P7_CaData00005',
    'SamplePupTrace_14',
    '111617_P8_CaData_00100005',
    '111617_P8_CaData_00100006',
    '121117_Animal_2_P8_CaData00004',
    '121117_Animal_2_P8_CaData00005'
};

masks = {
    '051118_P1_Ani2_Ca_00002-dFof_mask.mat',
    '120817_animal1_P1_CaData00001-dFof_mask.mat',
    '050218_P2_Ani2_Ca_00001-dFof_mask.mat',
    '050218_P2_Ani2_Ca_00002-dFof_mask.mat',
    '042318_P3_Ani1_00002-dFof_mask.mat',
    '042318_P3_Ani1_00004-dFof_mask.mat',
    '042318_P3_Ani1_00005-dFof_mask.mat',
    '042318_P3_Ani1_00007-dFof_mask.mat',
    '042518_P5_Ani1_00001-dFof_mask.mat',
    '042518_P5_Ani3_Ca_00001_67_98281-dFof_mask.mat',
    '120617_animal4_P7_CaData00019-dFof_mask.mat',
    '120617_animal4_P7_CaData00020-dFof_mask.mat',
    '121017_Animal_3_P7_CaData00001-dFof_mask.mat',
    '121017_Animal_3_P7_CaData00002-dFof_mask.mat',
    '121017_Animal_3_P7_CaData00005-dFof_mask.mat',
    'SamplePupTrace_14_mask.mat',
    '111617_P8_CaData_00100005-dFof_mask.mat',
    '111617_P8_CaData_00100006-dFof_mask.mat',
    '121117_Animal_2_P8_CaData00004-dFof_mask.mat',
    '121117_Animal_2_P8_CaData00005-dFof_mask.mat' 
};
    
basedir = '../data/figure7_ftle_data/';
maskdir = '../data/figure7_ftle_data/masks/';
savedir = './figure7/'; mkdir(savedir)

figsize = [0 0 240 180];
forward_color = [255 165 0]/255; backward_color = [134,0,212]/255;


% To exactly recreate the images in the manuscript please download and add ColorBrewer and Colorcet to the path:
%   https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
%   https://peterkovesi.com/projects/colourmaps/ 
%   AND uncomment below
% cmap = colorcet('L1');

cmap = colormap('gray');

for dir_i = 1:numel(dir_list)
    curr_dir = strcat(basedir,'ftle/', dir_list{dir_i});
    files = ls(curr_dir);
    files = strsplit(files,'.mat');
    load(strcat(maskdir,masks{dir_i}));

    if exist('mask','var')
        Mask = mask;
    end
    Msk=discard_zero_cols(double(Mask));
    Msk=BinStack(Msk,3);
    Msk = logical(Msk);
    

    mkdir(strcat(savedir,'flow/', dir_list{dir_i},'/'));
    mkdir(strcat(savedir,'quantScore/', dir_list{dir_i},'/'));

    for file_i = 1:numel(files)
        parts = strsplit(strtrim(files{file_i}),'_');

        if numel(parts{1}>0)
            curr_file = strcat(curr_dir, '/', strtrim(files{file_i}),'.mat');
            load(curr_file);
            
            if numel(parts{end}) == numel('sleep')

                if sum(size(mean_ftle.f))
                    % compute standard 0.93 thresholding, save FLOW portrait 
                    [for_FLOW, back_FLOW] = createFLOWPortrait(mean_ftle.f, mean_ftle.r, 0.93);

                    % % load mean_dFoF
                    load(strcat(basedir,'mean_dFoF/mean_dFoF',dir_list{dir_i},'-dFof.mat'));
                    % check dimensions of mean dFoF with flow
                    if ~isequal(size(for_FLOW),size(mean_dFoF)) % remove cols
                        mean_dFoF = mean_dFoF .* Mask;
                        mean_dFoF = discard_zero_cols(mean_dFoF);
                        mean_dFoF_ = BinStack(mean_dFoF,3)
                    end

                    if ~isequal(size(Msk), size(mean_dFoF))
                        1
                        Msk = BinStack(double(Mask),3);
                        Msk = discard_zero_cols(logical(Msk));
                        Msk = logical(Msk);
                    end

                    filename=strcat(savedir,'flow/', dir_list{dir_i},'/', strtrim(files{file_i}),'FLOW.png');
                    saveFlowPortrait(for_FLOW, back_FLOW, mean_dFoF, filename, true, Msk, forward_color, backward_color, cmap, figsize)
                    
                    % calculate and save ridgeCountScore for consolidation analysis
                    output = zeros(3,1);
                    [output(1), output(2), output(3)] = RidgeCountScore(for_FLOW, back_FLOW);
                    filename=strcat(savedir,'quantScore/', dir_list{dir_i},'/', strtrim(files{file_i}),'score.mat');
                    save(filename, 'output');
                end     
            end
        end
    end
    clear mask Mask
end

close all

%%%% Plot ridge count score and do statistical analysis %%%%
datapath = [savedir, 'quantScore/'];

% get cell array of all files
files = dir(datapath);
recordingFolders = {files([files.isdir]).name};

ages = [];
recordingID = {};
vals = [];

recording_mean93 = [];
recording_age = [];

% loop over all recording folders
for recIdx = 1:numel(recordingFolders);
    % get all file names without .
    currentRecording = dir(strcat(datapath, '/', recordingFolders{recIdx}));
    quantFiles = {currentRecording.name};
    quantFiles = quantFiles(~ismember(quantFiles,{'.','..'}));

    recording_scores93 = [];
    
    % loop over files in the current directory
    for fileIdx = 1:numel(quantFiles)
       currID = quantFiles{fileIdx};
       % only analyze sleep data
        if strfind(currID, 'sleep')
            recordingID{end+1} = currID;
            currData = load(strcat(datapath, '/', recordingFolders{recIdx}, ...
                '/', currID));
            
            % store age
            currAgeIdx = strfind(quantFiles{fileIdx},'P');
            if strfind(quantFiles{fileIdx}, 'SamplePupTrace_14')
                    age = 7
            else 
                    age = str2num(currID(currAgeIdx+1)); % only works since all ages are single digit
            end
            ages =  [ages age];

            % store age once per animal
            if fileIdx ==  1
                recording_age = [recording_age age];
            end
   
            vals = [vals currData.output];
            recording_scores93 = [recording_scores93 currData.output];
        end
    end 

    % store mean score for each recording
    recording_mean93 = [recording_mean93 mean(recording_scores93,2)];
end

% sort all data by age
[sortedAge, permutation] = sort(ages);
sortedVals = vals(:,permutation);
ages = unique(ages);

% compute mean and std of ages
groupedVals = cell(1,numel(ages));
stats = zeros(2,numel(ages),3);
for i = 1:numel(ages)
    currData = sortedVals(:,sortedAge==ages(i))
    groupedVals{i} = currData
    meanTemp = mean(currData,2,'omitnan');
    stdTemp = std(currData,0,2,'omitnan')/sqrt(size(currData, 2));
    for j =1:3
        stats(1,i,j) = meanTemp(j);
        stats(2,i,j) = stdTemp(j);
    end
end

% Plot %
msize=25; marker='_'; linewid = 2;

% 
figure('Renderer','painters','Position',[0 0 400 500])
scatter(sortedAge, sortedVals(1,:),msize,marker); hold on
errorbar(ages, stats(1,:,1), stats(2,:,1), 'k.','MarkerSize',25,'LineWidth',linewid)
yticks([0, 0.005, 0.01, 0.015, 0.02, 0.025]); 
xticks([1,2,3,5,7,8]); xticklabels({'P1','P2','P3','P5','P7','P8'});
xlim([0, 9]); ylim([0 0.03]); hold off
exportgraphics(gcf, [savedir, 'ridgeCountScore_forward.pdf'])

figure('Renderer','painters','Position',[0 0 400 500])
scatter(sortedAge, sortedVals(2,:),msize,marker); hold on
errorbar(ages, stats(1,:,2), stats(2,:,2),'k.','MarkerSize',25,'LineWidth',linewid)
yticks([0.0015, 0.003, 0.0045, 0.006, 0.0075, 0.009]);
xticks([1,2,3,5,7,8]); xticklabels({'P1','P2','P3','P5','P7','P8'});
xlim([0, 9]); ylim([0 0.0085]);  hold off
exportgraphics(gcf, [savedir, 'ridgeCountScore_backward.pdf'])

figure('Renderer','painters','Position',[0 0 400 500])
scatter(sortedAge, sortedVals(3,:),msize,marker); hold on
errorbar(ages, stats(1,:,3), stats(2,:,3), 'k.','MarkerSize',25,'LineWidth',linewid)
yticks([0.0015, 0.003, 0.0045, 0.006, 0.0075, 0.009]); 
xticks([1,2,3,5,7,8]); xticklabels({'P1','P2','P3','P5','P7','P8'})
xlim([0, 9]); ylim([0 0.008]); hold off
exportgraphics(gcf, [savedir, 'ridgeCountScore_combined.pdf'])

close all


% Perform statistical test of mean of P1-P3 and P5-P8
% Each column correponds to forward, backward, combined score respectively
young = sortedVals(:,sortedAge <= 3)';
old = sortedVals(:,sortedAge > 3)';
% compute paired t-test
% does test along the cols of young and old 
[h, p, ci, stats] = ttest2(young, old); 

% print to screen
formatSpecFor = 'The p-value for forward is %10.9f\n';
formatSpecBack = 'The p-value for backward is %10.9f\n';
formatSpecCombined = 'The p-value for combined is %10.9f\n';

fid = fopen([savedir, 'results.txt'],'w');
fprintf(fid, formatSpecFor, p(1));
fprintf(fid, formatSpecBack, p(2));
fprintf(fid, formatSpecCombined, p(3));
fclose(fid);

%%%% FUNCTIONS %%%%
function [modified_data, zero_idxs] = discard_zero_cols(data,discard_rows)
	% convert matrix to a logical array and take the sum of a row
	% a full row of zeros will come out as 0
	if nargin < 2
		discard_rows=false;
	end

	% discard cols
	zero_cols = sum(logical(data),1); 
	always_zero_cols = sum(zero_cols,3);  % add up in time to ensure only discarding constanly zero cols
	modified_data = data(:,find(always_zero_cols),:);  % extract non-zero columns and return

	% discard rows
	if discard_rows
		zero_rows = sum(logical(data),2);
		always_zero_rows = sum(zero_rows,3);  % add up in time to ensure only discarding constanly zero cols
		modified_data = modified_data(find(always_zero_rows),:,:);  % extract non-zero columns and return
	end
end

function [Binned] = BinStack(Data, BinSize)
    % Author: Dennis Tabuenna (UW)
    % 1 = no binning
    % disp('Binning...')
    
    
    Ystop = ceil(size(Data,1)/BinSize)*BinSize;
    Xstop = ceil(size(Data,2)/BinSize)*BinSize;
    
    Yrange = (size(Data,1)+1):Ystop;
    XRange = (size(Data,2)+1):Xstop;
    
    Data(Yrange, :, :) = NaN;
    Data(:, XRange, :) = NaN;
    
    for i = 1:BinSize
       Temp = Data(i:BinSize:Ystop, i:BinSize:Xstop, :, :);
       SubArrays(1:size(Temp,1),1:size(Temp,2),1:size(Temp,3),i) = Temp;
    end
    
    Binned = nanmean(SubArrays, 4);
    % disp('     done.')
end

function [forward, backward, forOrBack] = RidgeCountScore(for_FLOW, back_FLOW)
    [forw, back, fOrb] = localRidgeCount(for_FLOW, back_FLOW);
    [forArea, backArea, fOrbArea] = meanRidgeArea(for_FLOW, back_FLOW);

    forward = forw / forArea;
    backward = back / backArea;
    forOrBack = fOrb / fOrbArea;
end

function [forward, backward, forOrBack] = localRidgeCount(for_FLOW, back_FLOW)
    [~,numberForward] = bwlabel(for_FLOW);
    [~,numberBackward] = bwlabel(back_FLOW);
    [~,numberForOrBack] = bwlabel(for_FLOW | back_FLOW);
    
    forward = numberForward /  sum(sum(for_FLOW));
    backward = numberBackward / sum(sum(back_FLOW));
    forOrBack = numberForOrBack /  sum(sum(for_FLOW | back_FLOW));
    
    if isnan(forward)
        forward = 0;
    end
    
    if isnan(backward)
        backward = 0;
    end
    
    if isnan(forOrBack)
        forOrBack = 0;
    end
end

function [forward, backward, forOrBack] = meanRidgeArea(for_FLOW, back_FLOW)    
    [forInfo, IM_f] = bwferet(for_FLOW, 'MaxFeretProperties');
    [backInfo, IM_b] = bwferet(back_FLOW, 'MaxFeretProperties');
    [forOrBackInfo, IM_fORb] = bwferet((for_FLOW|back_FLOW), 'MaxFeretProperties');
    
    forward = mean(forInfo.MaxDiameter,'omitnan');
    backward = mean(backInfo.MaxDiameter);
    forOrBack = mean(forOrBackInfo.MaxDiameter);
end  