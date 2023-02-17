%% NIRS_ANALYSIS
% Loads the nirs data, preprocesses the data and performs a GLM analysis.
% The data is saved for later statistical analysis (nirs_statistics.m) or 
% visualization (nirs_plot.m).
%
% INPUT:
% - *.snirf files (bids_raw: nirs)
% - *events.tsv files (bids_raw: nirs)
% - TIMING.mat (scripts: nirs)
%
% OUTPUT:
% - gamma.m: output of the GLM analysis with gamma function for each participant (for statistics)
% - gaussian.m: output of the GLM analysis with gaussians for each participant (for visualization)
% 
% DEPENDENCIES:
% - Homer3 toolbox
% - custom_Homerfunct: custom Homer 3 functions

%% Load required packages & custom functions
cd     C:\Users\helen\Documents\MATLAB\matlab_toolboxes\Homer3
setpaths
rmpath(genpath('C:\Users\helen\Documents\MATLAB\matlab_toolboxes\fieldtrip')) % conflict between homer and fieldtrip toolbox
clc
addpath('C:\Users\helen\Documents\1.Projects\automaticity_fnirs\data\scripts\final\nirs\custom_Homerfct')

load('C:\Users\helen\Documents\1.Projects\automaticity_fnirs\data\scripts\final\nirs\TIMING.mat') % contains the relative onset of the finger presses/foot lifting (compared to the red cross) and duration of the trials (extracted from behavioral data)

%% Load data
raw_dir = 'F:\automaticity_fnirs\bids_raw\'; % data directory
proc_dir = 'C:\Users\helen\Documents\1.Projects\automaticity_fnirs\data\processed\final\nirs'; % directory to save processed data
sub_num=[6 19 21 27 31 33 34 38 43 44 52 53 63 67 69 74 75 77 84 87 88 90 95];

tasks = {'fingerauto', 'fingernonauto', 'footauto', 'footnonauto'};
tasks_timing = {'handauto', 'handnonauto', 'footauto', 'footnonauto'}; % in TIMING hand is used instead of finger

clear data_nirs dc_runs mlActRuns stimRuns

%% Process nirs data
for s = 1:length(sub_num)
  fprintf('\n===== Processing sub-%02d =====\n', sub_num(s))
for t = 1:numel(tasks)
  fprintf('--- task %s ---\n', tasks{t})
  % load the data
  filename = fullfile(raw_dir, sprintf('sub-%02d', sub_num(s)), 'nirs', sprintf('sub-%02d_task-%s_nirs.snirf', sub_num(s), tasks{t}));
  data_nirs = SnirfClass(filename); 
  
  % get optical densities
  dod = data_nirs.data.copy; % artinis fNIRS data is exported as optical densities
  
  % downsampling to 5 Hz
  dod_rs = dod;
  dod_rs.dataTimeSeries = downsample(dod.dataTimeSeries, 10); % original sampling rate is 50 Hz 
  nsmp = size(dod_rs.dataTimeSeries,1);
  dod_rs.time = transpose((0:(nsmp-1))/5);
  
  % remove channels with SCI <0.75
  probe = data_nirs.probe;
  mlActMan = cell(1,1);
  tIncMan = cell(1,1);
  threshold = 0.75;
  mlActRuns{t} = hmrR_PruneChannelsSCI_HC(dod_rs, probe, mlActMan, tIncMan, threshold); % custom function
  
  % load events of this run & exclude test trials & erroneous trials
  events_file = fullfile(raw_dir, sprintf('sub-%02d', sub_num(s)), 'nirs', sprintf('sub-%02d_task-%s_events.tsv', sub_num(s), tasks{t}));
  events = readtable(events_file, 'FileType', 'text', 'Delimiter', 'tab', 'ReadVariableNames', 1, 'TreatAsEmpty', 'n/a');
  stimonset_idx = find(strcmp(events.shown_stimulus, 'red X')); % sign to start tapping
  timing = TIMING(sub_num(s)).(tasks_timing{t}); % TIMING contains the relative offset of the first finger tap/first lifing of the leg compared to the appearance of the red X
  try
    taskonset = events.onset(stimonset_idx) + timing.onset;
  catch
    taskonset = events.onset(stimonset_idx) + 0.65*ones(11,1); % use the mean offset instead (e.g., for sub-67 motion data of the leg was missing)
  end
  excl = events.included_trial(stimonset_idx)==0; % trials to exclude
  
  % convert events to snirf format
  clear stim
  stim = StimClass();
  for i=1:length(tasks)
    stim(i).name = tasks{i};
  end
  stim(t).data = [taskonset 7*ones(size(taskonset)) ones(size(taskonset))];
  stim(t).states = [taskonset ones(size(taskonset))];
  stim(t).states(excl,2) = -2; 
  stimRuns{t} = stim; % concatenate stim of all runs
  
  % low pass filter: 0.25 Hz, 3rd order Butterworth
  dod_filt = hmrR_BandpassFilt(dod_rs, 0, 0.25);
  
  % short channel regression
  tIncAuto = cell(1,1);
  dod_regr = hmrR_ShortChannelRegression_HC(dod_filt, probe, mlActRuns{t}, tIncAuto, 15, 0, 2); % custom function
  
  % conversion to hemoglobin concentration changes
  dc = hmrR_OD2Conc(dod_regr, probe, [1 1]);
  dcRuns{t} = dc; % concatenate dc of all runs
end

%% GLM analysis
fprintf('--- GLM analysis ---\n')
for j=1:2
  % first GLM analysis for each condition separately (j=1), then GLM analysis for
  % finger and foot tasks merged (j=2)
  if j==1 % auto vs non-auto
    save_dir = fullfile(proc_dir, sprintf('auto_vs_nonauto'), sprintf('sub-%02d', sub_num(s)));
    mkdir(save_dir);
  elseif j==2 % finger vs foot
    save_dir = fullfile(proc_dir, sprintf('finger_vs_foot'), sprintf('sub-%02d', sub_num(s)));
    mkdir(save_dir);
    % merge the finger and foot events
    for t = 1:numel(tasks)
      clear stim
      stim = StimClass();
      stim(1).name = 'finger';
      stim(2).name = 'foot'; 
      if contains(tasks(t), 'finger')
        stim(1).data = stimRuns{t}(t).data;
        stim(1).states = stimRuns{t}(t).states;
      elseif contains(tasks(t), 'foot')
        stim(2).data = stimRuns{t}(t).data;
        stim(2).states = stimRuns{t}(t).states;
      end
      stimRuns{t} = stim; 
    end
  end
  
  % GLM analysis with gamma function(for each participant)
  AauxRuns = cell(1,4);
  tIncAutoRuns = cell(1,4);
  rcMapRuns = cell(1,4);
  trange = [-10 20];
  glmSolveMethod = 2;
  idxBasis = 2; 
  paramsBasis = [0.1 3 8 1.8 3 8];
  rhoSD_ssThresh = -1; % omit short channel regression because already performed earlier (because of multicollinearity)
  % note that sub-69 had a short channel with inaccurate optode positions,
  % and was estimated with a distance of 0 (so rhoSD_ssThresh = 0 would not work for
  % this subject)
  flagNuisanceRMethod = 0;
  driftOrder = 3;
  c_vector = [1 -1 0 0];
%   [yavg, yavgstd, nTrials, ynew, yresid, ysum2, beta, R, hmrstats] = hmrS_GLM_HC(dcRuns, stimRuns, probe, mlActRuns, AauxRuns, tIncAutoRuns, rcMapRuns, trange, glmSolveMethod, idxBasis, paramsBasis, rhoSD_ssThresh, flagNuisanceRMethod, driftOrder, c_vector);
%   % save beta for each subject
%   gamma.beta = beta;
%   gamma.hmrstats = hmrstats;
%   gamma.nTrials = nTrials;
%   filename = fullfile(save_dir, 'gamma');
%   save(filename, 'gamma');
  
  % GLM analysis with gaussians(for each participant) (for visualization)
  idxBasis = 1;
  paramsBasis = [1 1 0 0 0 0]; % set to default for gaussians
  [yavg, yavgstd, nTrials, ynew, yresid, ysum2, beta, R, hmrstats] = hmrS_GLM_HC(dcRuns, stimRuns, probe, mlActRuns, AauxRuns, tIncAutoRuns, rcMapRuns, trange, glmSolveMethod, idxBasis, paramsBasis, rhoSD_ssThresh, flagNuisanceRMethod, driftOrder, c_vector);
  gaussian.yavg = yavg;
  gaussian.hmrstats = hmrstats;
  gaussian.nTrials = nTrials;
  filename = fullfile(save_dir, 'gaussian');
  save(filename, 'gaussian');
end
end


