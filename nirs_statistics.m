%% NIRS_STATISTICS
% performs the t-test statistics based on the GLM output with gamma
% functions. 
%
% INPUT:
% - gamma.mat (processed: nirs)
% 
% OUTPUT:
% - stat_*.tsv: output table with the statistics for each condition and
% chromophore
% - N_by_channel.mat: overview of the number of participants for each
% channel

%% FINGER vs LEG
proc_dir = 'C:\Users\helen\Documents\1.Projects\automaticity_fnirs\data\processed\final\nirs\finger_vs_foot';
cd(proc_dir)
sub = dir('sub-*');

% create empty matrices
BETA_1 = nan(length(sub), 36, 2); % [sub x chan x conc]
BETA_2 = nan(length(sub), 36, 2);

% load beta values for each participant
for s=1:length(sub)
  load([sub(s).name '/gamma.mat'])
  beta_1 = permute(gamma.beta{1}(1, :, :, 1), [1 3 2]); % [#coefficients x HbX x #Channels x #conditions]
  BETA_1(s, :, :) = beta_1;
  beta_2 = permute(gamma.beta{1}(1, :, :, 2), [1 3 2]);
  BETA_2(s, :, :) = beta_2;
end

% replace zeros by nans (Homer attribute zero's to nan values)
BETA_1(find(BETA_1==0))=nan;
BETA_2(find(BETA_2==0))=nan;

% count number of participants for each channel
N_by_channel = sum(~isnan(BETA_1(:,:,2)));
if N_by_channel ~= sum(~isnan(BETA_2(:,:,2)))
  error
end
save('N_by_channel')

% channels in the order that they appear in the manuscript
chan = [5 4 12 11 8 13 9]; % order of channels in manuscript
chan_names = [1:7]'; % names according to manuscript

% FINGER VS BASELINE
[h,p,ci,stats] = ttest(BETA_1(:,chan,:));
p_fdr{1} = mafdr(permute(p(1,:,1),[2, 3,1]), 'BHFDR', true); % FDR correction of p-values
p_fdr{2} = mafdr(permute(p(1,:,2),[2,3,1]), 'BHFDR', true);
d{1} = permute(stats.tstat(1, :, 1), [2 3 1])./sqrt(N_by_channel(chan))'; % cohen's d
d{2} = permute(stats.tstat(1, :, 2), [2 3 1])./sqrt(N_by_channel(chan))';
chrom = {'HbO', 'HbR'};
for c=1:2 % HbO vs HbR
  mean = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanmean(BETA_1(:,chan,c))');
  sd = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanstd(BETA_1(:,chan,c))');
  CI = arrayfun(@(x,y) {sprintf('[%.1f %.1f]', x, y)}, 10^6*ci(1,:,c)', 10^6*ci(2,:,c)');
  t_value = arrayfun(@(x) {sprintf('%.1f', x)}, stats.tstat(:,:,c)');
  p_value = arrayfun(@(x) {sprintf('%.3f', x)}, p_fdr{c});
  p_value(p_fdr{c}<0.001) = {'<0.001'};
  Cohens_d = arrayfun(@(x) {sprintf('%.1f', x)}, d{c});
t=table(chan_names, mean, sd, CI,...
  t_value, stats.df(:,:,c)', p_value, Cohens_d,...
  'VariableNames',{'Channel', 'M', 'SD', 'CI', 't', 'df', 'p', 'd'})
tablefile = fullfile(proc_dir, sprintf('stat_finger_%s.tsv', chrom{c}));
writetable(t, tablefile, 'FileType', 'text','Delimiter', 'tab')
end

% LEG VS BASELINE
[h,p,ci,stats] = ttest(BETA_2(:,chan,:));
p_fdr{1} = mafdr(permute(p(1,:,1),[2, 3,1]), 'BHFDR', true); % FDR correction of p-values
p_fdr{2} = mafdr(permute(p(1,:,2),[2,3,1]), 'BHFDR', true);
d{1} = permute(stats.tstat(1, :, 1), [2 3 1])./sqrt(N_by_channel(chan))'; % cohen's d
d{2} = permute(stats.tstat(1, :, 2), [2 3 1])./sqrt(N_by_channel(chan))';
for c=1:2 % HbO vs HbR
  mean = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanmean(BETA_2(:,chan,c))');
  sd = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanstd(BETA_2(:,chan,c))');
  CI = arrayfun(@(x,y) {sprintf('[%.1f %.1f]', x, y)}, 10^6*ci(1,:,c)', 10^6*ci(2,:,c)');
  t_value = arrayfun(@(x) {sprintf('%.1f', x)}, stats.tstat(:,:,c)');
  p_value = arrayfun(@(x) {sprintf('%.3f', x)}, p_fdr{c});
  p_value(p_fdr{c}<0.001) = {'<0.001'};
  Cohens_d = arrayfun(@(x) {sprintf('%.1f', x)}, d{c});
t=table(chan_names, mean, sd, CI,...
  t_value, stats.df(:,:,c)', p_value, Cohens_d,...
  'VariableNames',{'Channel', 'M', 'SD', 'CI', 't', 'df', 'p', 'd'})
tablefile = fullfile(proc_dir, sprintf('stat_leg_%s.tsv', chrom{c}));
writetable(t, tablefile, 'FileType', 'text','Delimiter', 'tab')
end

% FINGER VS LEG
[h,p,ci,stats] = ttest(BETA_1(:,chan,:), BETA_2(:,chan,:));
p_fdr{1} = mafdr(permute(p(1,:,1),[2, 3,1]), 'BHFDR', true); % FDR correction of p-values
p_fdr{2} = mafdr(permute(p(1,:,2),[2,3,1]), 'BHFDR', true);
d{1} = permute(stats.tstat(1, :, 1), [2 3 1])./sqrt(N_by_channel(chan))'; % cohen's d
d{2} = permute(stats.tstat(1, :, 2), [2 3 1])./sqrt(N_by_channel(chan))';
for c=1:2 % HbO vs HbR
  mean = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*(nanmean(BETA_1(:,chan,c))-nanmean(BETA_2(:,chan,c)))');
  sd = arrayfun(@(x) {sprintf('%.1f', x)}, (10^6*(nanstd(BETA_1(:,chan,c))+nanstd(BETA_2(:,chan,c)))/2)');
  CI = arrayfun(@(x,y) {sprintf('[%.1f %.1f]', x, y)}, 10^6*ci(1,:,c)', 10^6*ci(2,:,c)');
  t_value = arrayfun(@(x) {sprintf('%.1f', x)}, stats.tstat(:,:,c)');
  p_value = arrayfun(@(x) {sprintf('%.3f', x)}, p_fdr{c});
  p_value(p_fdr{c}<0.001) = {'<0.001'};
  Cohens_d = arrayfun(@(x) {sprintf('%.1f', x)}, d{c});
t=table(chan_names, mean, sd, CI,...
  t_value, stats.df(:,:,c)', p_value, Cohens_d,...
  'VariableNames',{'Channel', 'M', 'SD', 'CI', 't', 'df', 'p', 'd'})
tablefile = fullfile(proc_dir, sprintf('stat_fingervsleg_%s.tsv', chrom{c}));
writetable(t, tablefile, 'FileType', 'text', 'Delimiter', 'tab')
end

%% AUTO vs NONAUTO
proc_dir = 'C:\Users\helen\Documents\1.Projects\automaticity_fnirs\data\processed\final\nirs\auto_vs_nonauto';
cd(proc_dir)
sub = dir('sub-*');

% create empty matrices
BETA_1 = nan(length(sub), 36, 2); % [sub x chan x conc]
BETA_2 = nan(length(sub), 36, 2);
BETA_3 = nan(length(sub), 36, 2);
BETA_4 = nan(length(sub), 36, 2);

% load beta values for each participant
sub = dir('sub-*');
for s=1:length(sub)
  load([sub(s).name '/gamma.mat'])
  beta_1 = permute(gamma.beta{1}(1, :, :, 1), [1 3 2]); % [#coefficients x HbX x #Channels x #conditions]
  BETA_1(s, :, :) = beta_1;
  beta_2 = permute(gamma.beta{1}(1, :, :, 2), [1 3 2]);
  BETA_2(s, :, :) = beta_2;
  beta_3 = permute(gamma.beta{1}(1, :, :, 3), [1 3 2]);
  BETA_3(s, :, :) = beta_3;
  beta_4 = permute(gamma.beta{1}(1, :, :, 4), [1 3 2]);
  BETA_4(s, :, :) = beta_4;
end

% replace zeros by nans (Homer attribute zero's to nan values)
BETA_1(find(BETA_1==0))=nan;
BETA_2(find(BETA_2==0))=nan;
BETA_3(find(BETA_3==0))=nan;
BETA_4(find(BETA_4==0))=nan;

% count number of participants for each channel
N_by_channel = sum(~isnan(BETA_1(:,:,2)));
if N_by_channel ~= sum(~isnan(BETA_2(:,:,2))) | N_by_channel ~= sum(~isnan(BETA_3(:,:,2))) | N_by_channel ~= sum(~isnan(BETA_4(:,:,2)))
  error
end
save('N_by_channel')

% channels in the order that they appear in the manuscript
chan_manuscript = [6 2; 4 2; 6 4; 4 4; 4 3; 7 4; 7 3; 4 1; 2 1; 2 3; 12 6; 12 8; 14 8; 14 6;...
  17 9; 17 10; 18 10; 18 9; 10 5; 10 7; 15 7; 19 11; 19 12; 20 12]; % channels (sources-detectors) in the order that they are presented in the manuscript
[~, chan] = ismember(chan_manuscript,gamma.hmrstats.ml, 'rows');
chan_names = [1:24]'; % names according to manuscript

% FINGER: AUTO VS BASELINE
[h,p,ci,stats] = ttest(BETA_1(:,chan,:));
p_fdr{1} = mafdr(permute(p(1,:,1),[2, 3,1]), 'BHFDR', true); % FDR correction of p-values
p_fdr{2} = mafdr(permute(p(1,:,2),[2,3,1]), 'BHFDR', true);
d{1} = permute(stats.tstat(1, :, 1), [2 3 1])./sqrt(N_by_channel(chan))'; % cohen's d
d{2} = permute(stats.tstat(1, :, 2), [2 3 1])./sqrt(N_by_channel(chan))';
chrom = {'HbO', 'HbR'};
for c=1:2 % HbO vs HbR
  mean = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanmean(BETA_1(:,chan,c))');
  sd = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanstd(BETA_1(:,chan,c))');
  CI = arrayfun(@(x,y) {sprintf('[%.1f %.1f]', x, y)}, 10^6*ci(1,:,c)', 10^6*ci(2,:,c)');
  t_value = arrayfun(@(x) {sprintf('%.1f', x)}, stats.tstat(:,:,c)');
  p_value = arrayfun(@(x) {sprintf('%.3f', x)}, p_fdr{c});
  p_value(p_fdr{c}<0.001) = {'<0.001'};
  Cohens_d = arrayfun(@(x) {sprintf('%.1f', x)}, d{c});
t=table(chan_names, mean, sd, CI,...
  t_value, stats.df(:,:,c)', p_value, Cohens_d,...
  'VariableNames',{'Channel', 'M', 'SD', 'CI', 't', 'df', 'p', 'd'})
tablefile = fullfile(proc_dir, sprintf('stat_auto_finger_%s.tsv', chrom{c}));
writetable(t, tablefile, 'FileType', 'text','Delimiter', 'tab')
end

% FINGER: NONAUTO VS BASELINE
[h,p,ci,stats] = ttest(BETA_2(:,chan,:));
p_fdr{1} = mafdr(permute(p(1,:,1),[2, 3,1]), 'BHFDR', true); % FDR correction of p-values
p_fdr{2} = mafdr(permute(p(1,:,2),[2,3,1]), 'BHFDR', true);
d{1} = permute(stats.tstat(1, :, 1), [2 3 1])./sqrt(N_by_channel(chan))'; % cohen's d
d{2} = permute(stats.tstat(1, :, 2), [2 3 1])./sqrt(N_by_channel(chan))';
chrom = {'HbO', 'HbR'};
for c=1:2 % HbO vs HbR
  mean = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanmean(BETA_2(:,chan,c))');
  sd = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanstd(BETA_2(:,chan,c))');
  CI = arrayfun(@(x,y) {sprintf('[%.1f %.1f]', x, y)}, 10^6*ci(1,:,c)', 10^6*ci(2,:,c)');
  t_value = arrayfun(@(x) {sprintf('%.1f', x)}, stats.tstat(:,:,c)');
  p_value = arrayfun(@(x) {sprintf('%.3f', x)}, p_fdr{c});
  p_value(p_fdr{c}<0.001) = {'<0.001'};
  Cohens_d = arrayfun(@(x) {sprintf('%.1f', x)}, d{c});
t=table(chan_names, mean, sd, CI,...
  t_value, stats.df(:,:,c)', p_value, Cohens_d,...
  'VariableNames',{'Channel', 'M', 'SD', 'CI', 't', 'df', 'p', 'd'})
tablefile = fullfile(proc_dir, sprintf('stat_nonauto_finger_%s.tsv', chrom{c}));
writetable(t, tablefile, 'FileType', 'text','Delimiter', 'tab')
end

% FINGER: AUTO VS NONAUTO
[h,p,ci,stats] = ttest(BETA_1(:,chan,:), BETA_2(:,chan,:));
p_fdr{1} = mafdr(permute(p(1,:,1),[2, 3,1]), 'BHFDR', true); % FDR correction of p-values
p_fdr{2} = mafdr(permute(p(1,:,2),[2,3,1]), 'BHFDR', true);
d{1} = permute(stats.tstat(1, :, 1), [2 3 1])./sqrt(N_by_channel(chan))'; % cohen's d
d{2} = permute(stats.tstat(1, :, 2), [2 3 1])./sqrt(N_by_channel(chan))';
for c=1:2 % HbO vs HbR
  mean = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*(nanmean(BETA_1(:,chan,c))-nanmean(BETA_2(:,chan,c)))');
  sd = arrayfun(@(x) {sprintf('%.1f', x)}, (10^6*(nanstd(BETA_1(:,chan,c))+nanstd(BETA_2(:,chan,c)))/2)');
  CI = arrayfun(@(x,y) {sprintf('[%.1f %.1f]', x, y)}, 10^6*ci(1,:,c)', 10^6*ci(2,:,c)');
  t_value = arrayfun(@(x) {sprintf('%.1f', x)}, stats.tstat(:,:,c)');
  p_value = arrayfun(@(x) {sprintf('%.3f', x)}, p_fdr{c});
  p_value(p_fdr{c}<0.001) = {'<0.001'};
  Cohens_d = arrayfun(@(x) {sprintf('%.1f', x)}, d{c});
t=table(chan_names, mean, sd, CI,...
  t_value, stats.df(:,:,c)', p_value, Cohens_d,...
  'VariableNames',{'Channel', 'M', 'SD', 'CI', 't', 'df', 'p', 'd'})
tablefile = fullfile(proc_dir, sprintf('stat_autovsnonauto_finger_%s.tsv', chrom{c}));
writetable(t, tablefile, 'FileType', 'text','Delimiter', 'tab')
end

% LEG: AUTO VS BASELINE
[h,p,ci,stats] = ttest(BETA_3(:,chan,:));
p_fdr{1} = mafdr(permute(p(1,:,1),[2, 3,1]), 'BHFDR', true); % FDR correction of p-values
p_fdr{2} = mafdr(permute(p(1,:,2),[2,3,1]), 'BHFDR', true);
d{1} = permute(stats.tstat(1, :, 1), [2 3 1])./sqrt(N_by_channel(chan))'; % cohen's d
d{2} = permute(stats.tstat(1, :, 2), [2 3 1])./sqrt(N_by_channel(chan))';
chrom = {'HbO', 'HbR'};
for c=1:2 % HbO vs HbR
  mean = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanmean(BETA_3(:,chan,c))');
  sd = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanstd(BETA_3(:,chan,c))');
  CI = arrayfun(@(x,y) {sprintf('[%.1f %.1f]', x, y)}, 10^6*ci(1,:,c)', 10^6*ci(2,:,c)');
  t_value = arrayfun(@(x) {sprintf('%.1f', x)}, stats.tstat(:,:,c)');
  p_value = arrayfun(@(x) {sprintf('%.3f', x)}, p_fdr{c});
  p_value(p_fdr{c}<0.001) = {'<0.001'};
  Cohens_d = arrayfun(@(x) {sprintf('%.1f', x)}, d{c});
t=table(chan_names, mean, sd, CI,...
  t_value, stats.df(:,:,c)', p_value, Cohens_d,...
  'VariableNames',{'Channel', 'M', 'SD', 'CI', 't', 'df', 'p', 'd'})
tablefile = fullfile(proc_dir, sprintf('stat_auto_leg_%s.tsv', chrom{c}));
writetable(t, tablefile, 'FileType', 'text','Delimiter', 'tab')
end

% LEG: NONAUTO VS BASELINE
[h,p,ci,stats] = ttest(BETA_4(:,chan,:));
p_fdr{1} = mafdr(permute(p(1,:,1),[2, 3,1]), 'BHFDR', true); % FDR correction of p-values
p_fdr{2} = mafdr(permute(p(1,:,2),[2,3,1]), 'BHFDR', true);
d{1} = permute(stats.tstat(1, :, 1), [2 3 1])./sqrt(N_by_channel(chan))'; % cohen's d
d{2} = permute(stats.tstat(1, :, 2), [2 3 1])./sqrt(N_by_channel(chan))';
chrom = {'HbO', 'HbR'};
for c=1:2 % HbO vs HbR
  mean = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanmean(BETA_4(:,chan,c))');
  sd = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*nanstd(BETA_4(:,chan,c))');
  CI = arrayfun(@(x,y) {sprintf('[%.1f %.1f]', x, y)}, 10^6*ci(1,:,c)', 10^6*ci(2,:,c)');
  t_value = arrayfun(@(x) {sprintf('%.1f', x)}, stats.tstat(:,:,c)');
  p_value = arrayfun(@(x) {sprintf('%.3f', x)}, p_fdr{c});
  p_value(p_fdr{c}<0.001) = {'<0.001'};
  Cohens_d = arrayfun(@(x) {sprintf('%.1f', x)}, d{c});
t=table(chan_names, mean, sd, CI,...
  t_value, stats.df(:,:,c)', p_value, Cohens_d,...
  'VariableNames',{'Channel', 'M', 'SD', 'CI', 't', 'df', 'p', 'd'})
tablefile = fullfile(proc_dir, sprintf('stat_nonauto_leg_%s.tsv', chrom{c}));
writetable(t, tablefile, 'FileType', 'text','Delimiter', 'tab')
end

% LEG: AUTO VS NONAUTO
[h,p,ci,stats] = ttest(BETA_3(:,chan,:), BETA_4(:,chan,:));
p_fdr{1} = mafdr(permute(p(1,:,1),[2, 3,1]), 'BHFDR', true); % FDR correction of p-values
p_fdr{2} = mafdr(permute(p(1,:,2),[2,3,1]), 'BHFDR', true);
d{1} = permute(stats.tstat(1, :, 1), [2 3 1])./sqrt(N_by_channel(chan))'; % cohen's d
d{2} = permute(stats.tstat(1, :, 2), [2 3 1])./sqrt(N_by_channel(chan))';
for c=1:2 % HbO vs HbR
  mean = arrayfun(@(x) {sprintf('%.1f', x)}, 10^6*(nanmean(BETA_3(:,chan,c))-nanmean(BETA_4(:,chan,c)))');
  sd = arrayfun(@(x) {sprintf('%.1f', x)}, (10^6*(nanstd(BETA_3(:,chan,c))+nanstd(BETA_4(:,chan,c)))/2)');
  CI = arrayfun(@(x,y) {sprintf('[%.1f %.1f]', x, y)}, 10^6*ci(1,:,c)', 10^6*ci(2,:,c)');
  t_value = arrayfun(@(x) {sprintf('%.1f', x)}, stats.tstat(:,:,c)');
  p_value = arrayfun(@(x) {sprintf('%.3f', x)}, p_fdr{c});
  p_value(p_fdr{c}<0.001) = {'<0.001'};
  Cohens_d = arrayfun(@(x) {sprintf('%.1f', x)}, d{c});
t=table(chan_names, mean, sd, CI,...
  t_value, stats.df(:,:,c)', p_value, Cohens_d,...
  'VariableNames',{'Channel', 'M', 'SD', 'CI', 't', 'df', 'p', 'd'})
tablefile = fullfile(proc_dir, sprintf('stat_autovsnonauto_leg_%s.tsv', chrom{c}));
writetable(t, tablefile, 'FileType', 'text','Delimiter', 'tab')
end


