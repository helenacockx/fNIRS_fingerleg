%% NIRS_PLOT
% Plots the outcome of the GLM analysis with gaussians. Plots each channel
% separately with HbO in red and HbR in blue and a shaded area containing
% the 95% confidence intervals
%
% INPUT:
% gaussian.mat (processed: nirs)
% N_by_channel.mat (processed: nirs) (created by nirs_statistics.m)
%
% DEPENDENCIES:
% - Fieldtrip toolbox
% - (Homer3 toolbox)

%% Load Fieldtrip & Homer toolbox
cd     'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\MATLAB\matlab_toolboxes\Homer3'
setpaths

addpath     'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\MATLAB\matlab_toolboxes\fieldtrip'
ft_defaults;

%% FINGER vs LEG
proc_dir = 'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\automaticity_fnirs\data\processed\final\nirs\finger_vs_foot';
cd(proc_dir)
sub = dir('sub-*');

% create a new struct with the data of all participants
yAvgSubjs = cell(1, length(sub));
nTrialsSubjs = cell(1, length(sub));

% load yavg values for each participant
for s=1:length(sub)
  load([sub(s).name '/gaussian.mat'])
  yAvgSubjs{s} = gaussian.yavg;
  nTrialsSubjs{s} = gaussian.nTrials;
end

% average over all participants
[yAvg, nTrials] = hmrG_SubjAvg(yAvgSubjs, nTrialsSubjs);
[yAvgStd, yAvgStdErr] = hmrG_SubjAvgStd(yAvgSubjs);

N=length(sub); 
load('N_by_channel') % from nirs_statistics.m

ML=yAvg.measurementList;
num_con = 2;
for i = 1:3:length(yAvg.measurementList)/num_con
  n = N_by_channel((i+2)/3);
  t = yAvg.time;
  y1 = yAvg.dataTimeSeries(:, [i:i+1]);
  y1_std = yAvgStd.dataTimeSeries(:, [i:i+1]);
  y1_MOE = (1.96*y1_std)/sqrt(n); % z-score for 95% CI = 1.96
  y1_CI = [y1+y1_MOE y1-y1_MOE]; % [upperCI_O2Hb upperCI_HHb lowerCI_O2Hb lowerCI_HHb];
  % plot
  figure; 
  subplot(1, 2, 1); hold on; title('finger')
  ft_plot_box([0 8 -1.5*10^-5 3*10^-5], 'facecolor', [0.6 0.6 0.6], 'facealpha', 0.5, 'edgecolor', 'none')
  ft_plot_vector(t, y1_CI(:, [1 3])', 'highlight', ones(size(y1(:,1)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'r', 'facealpha', 0.3);
  ft_plot_vector(t, y1(:,1)', 'color', 'r', 'linewidth', 2)
  ft_plot_vector(t, y1_CI(:, [2 4])', 'highlight', ones(size(y1(:,2)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'b', 'facealpha', 0.3);
  ft_plot_vector(t, y1(:,2)', 'color', 'b', 'linewidth', 2)
  ylim([-0.7*10^-5 1.7*10^-5]);
  
  % find corresponding channel of the condition 2
  idx = find([ML(:).sourceIndex]==ML(i).sourceIndex & [ML(:).detectorIndex]==ML(i).detectorIndex & [ML(:).dataTypeIndex]==2);
  y2 = yAvg.dataTimeSeries(:, [idx(1) idx(2)]);
  y2_std = yAvgStd.dataTimeSeries(:, [idx(1) idx(2)]);
  y2_MOE = (1.96*y2_std)/sqrt(n); 
  y2_CI = [y2+y2_MOE y2-y2_MOE]; % [upperCI_O2Hb upperCI_HHb lowerCI_O2Hb lowerCI_HHb];
  % plot
   subplot(1,2,2); hold on; title('foot')
  ft_plot_box([0 8 -1.5*10^-5 3*10^-5], 'facecolor', [0.6 0.6 0.6], 'facealpha', 0.5, 'edgecolor', 'none')
  ft_plot_vector(t, y2_CI(:, [1 3])', 'highlight', ones(size(y2(:,1)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'r', 'facealpha', 0.3);
  ft_plot_vector(t, y2(:,1)', 'color', 'r', 'linewidth', 2)
  ft_plot_vector(t, y2_CI(:, [2 4])', 'highlight', ones(size(y2(:,2)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'b', 'facealpha', 0.3);
  ft_plot_vector(t, y2(:,2)', 'color', 'b', 'linewidth', 2)
  ylim([-0.7*10^-5 1.7*10^-5]);
  
  sgtitle(sprintf('S%d - D%d', ML(i).sourceIndex, ML(i).detectorIndex))
end

%% AUTO vs NON-AUTO
proc_dir = 'C:\Users\helen\OneDrive - Radboud Universiteit\Documenten\1.Projects\automaticity_fnirs\data\processed\final\nirs\auto_vs_nonauto';
cd(proc_dir)
sub = dir('sub-*');

% create a new struct with the data of all participants
yAvgSubjs = cell(1, length(sub));
nTrialsSubjs = cell(1, length(sub));

% load yavg values for each participant
for s=1:length(sub)
  load([sub(s).name '/gaussian.mat'])
  yAvgSubjs{s} = gaussian.yavg;
  nTrialsSubjs{s} = gaussian.nTrials;
end

% average over all participants
[yAvg, nTrials] = hmrG_SubjAvg(yAvgSubjs, nTrialsSubjs);
[yAvgStd, yAvgStdErr] = hmrG_SubjAvgStd(yAvgSubjs);

N=length(sub); 
load('N_by_channel') % from nirs_statistics.m

ML=yAvg.measurementList;
num_con = 4;
for i = 1:3:length(yAvg.measurementList)/num_con
  n = N_by_channel((i+2)/3);
  t = yAvg.time;
  y1 = yAvg.dataTimeSeries(:, [i:i+1]);
  y1_std = yAvgStd.dataTimeSeries(:, [i:i+1]);
  y1_MOE = (1.96*y1_std)/sqrt(n); % z-score for 95% CI = 1.96
  y1_CI = [y1+y1_MOE y1-y1_MOE]; % [upperCI_O2Hb upperCI_HHb lowerCI_O2Hb lowerCI_HHb];
  % plot
  figure; 
  subplot(2, 2, 1); hold on; title('finger auto')
  ft_plot_box([0 8 -1.5*10^-5 3*10^-5], 'facecolor', [0.6 0.6 0.6], 'facealpha', 0.5, 'edgecolor', 'none')
  ft_plot_vector(t, y1_CI(:, [1 3])', 'highlight', ones(size(y1(:,1)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'r', 'facealpha', 0.3);
  ft_plot_vector(t, y1(:,1)', 'color', 'r', 'linewidth', 2)
  ft_plot_vector(t, y1_CI(:, [2 4])', 'highlight', ones(size(y1(:,2)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'b', 'facealpha', 0.3);
  ft_plot_vector(t, y1(:,2)', 'color', 'b', 'linewidth', 2)
  ylim([-0.7*10^-5 1.7*10^-5]);
  
  % find corresponding channel of the condition 2
  idx = find([ML(:).sourceIndex]==ML(i).sourceIndex & [ML(:).detectorIndex]==ML(i).detectorIndex & [ML(:).dataTypeIndex]==2);
  y2 = yAvg.dataTimeSeries(:, [idx(1) idx(2)]);
  y2_std = yAvgStd.dataTimeSeries(:, [idx(1) idx(2)]);
  y2_MOE = (1.96*y2_std)/sqrt(n); 
  y2_CI = [y2+y2_MOE y2-y2_MOE]; % [upperCI_O2Hb upperCI_HHb lowerCI_O2Hb lowerCI_HHb];
  % plot
  subplot(2,2,3); hold on; title('finger nonauto')
  ft_plot_box([0 8 -1.5*10^-5 3*10^-5], 'facecolor', [0.6 0.6 0.6], 'facealpha', 0.5, 'edgecolor', 'none')
  ft_plot_vector(t, y2_CI(:, [1 3])', 'highlight', ones(size(y2(:,1)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'r', 'facealpha', 0.3);
  ft_plot_vector(t, y2(:,1)', 'color', 'r', 'linewidth', 2)
  ft_plot_vector(t, y2_CI(:, [2 4])', 'highlight', ones(size(y2(:,2)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'b', 'facealpha', 0.3);
  ft_plot_vector(t, y2(:,2)', 'color', 'b', 'linewidth', 2)
  ylim([-0.7*10^-5 1.7*10^-5]);
  
  % find corresponding channel of the condition 3
  idx = find([ML(:).sourceIndex]==ML(i).sourceIndex & [ML(:).detectorIndex]==ML(i).detectorIndex & [ML(:).dataTypeIndex]==3);
  y3 = yAvg.dataTimeSeries(:, [idx(1) idx(2)]);
  y3_std = yAvgStd.dataTimeSeries(:, [idx(1) idx(2)]);
  y3_MOE = (1.96*y3_std)/sqrt(n); 
  y3_CI = [y3+y3_MOE y3-y3_MOE]; % [upperCI_O2Hb upperCI_HHb lowerCI_O2Hb lowerCI_HHb];
  % plot
  subplot(2,2,2); hold on; title('foot auto')
  ft_plot_box([0 8 -1.5*10^-5 3*10^-5], 'facecolor', [0.6 0.6 0.6], 'facealpha', 0.5, 'edgecolor', 'none')
  ft_plot_vector(t, y3_CI(:, [1 3])', 'highlight', ones(size(y3(:,1)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'r', 'facealpha',0.3);
  ft_plot_vector(t, y3(:,1)', 'color', 'r', 'linewidth', 2)
  ft_plot_vector(t, y3_CI(:, [2 4])', 'highlight', ones(size(y3(:,2)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'b', 'facealpha', 0.3);
  ft_plot_vector(t, y3(:,2)', 'color', 'b', 'linewidth', 2)
  ylim([-1.5*10^-5 3*10^-5]);

  % find corresponding channel of the condition 4
  idx = find([ML(:).sourceIndex]==ML(i).sourceIndex & [ML(:).detectorIndex]==ML(i).detectorIndex & [ML(:).dataTypeIndex]==4);
  y4 = yAvg.dataTimeSeries(:, [idx(1) idx(2)]);
  y4_std = yAvgStd.dataTimeSeries(:, [idx(1) idx(2)]);
  y4_MOE = (1.96*y4_std)/sqrt(n); 
  y4_CI = [y4+y4_MOE y4-y4_MOE]; % [upperCI_O2Hb upperCI_HHb lowerCI_O2Hb lowerCI_HHb];
  % plot
  subplot(2,2,4); hold on; title('foot nonauto')
  ft_plot_box([0 8 -1.5*10^-5 3*10^-5], 'facecolor', [0.6 0.6 0.6], 'facealpha', 0.5, 'edgecolor', 'none')
  ft_plot_vector(t, y4_CI(:, [1 3])', 'highlight', ones(size(y4(:,1)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'r', 'facealpha', 0.3);
  ft_plot_vector(t, y4(:,1)', 'color', 'r', 'linewidth', 2)
  ft_plot_vector(t, y4_CI(:, [2 4])', 'highlight', ones(size(y4(:,2)')), 'highlightstyle', 'difference', 'color', 'none','facecolor', 'b', 'facealpha', 0.3);
  ft_plot_vector(t, y4(:,2)', 'color', 'b', 'linewidth', 2)
  ylim([-1.5*10^-5 3*10^-5]);
  
  sgtitle(sprintf('S%d - D%d', ML(i).sourceIndex, ML(i).detectorIndex))
end
