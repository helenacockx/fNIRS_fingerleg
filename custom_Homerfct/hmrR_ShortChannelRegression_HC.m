% SYNTAX:
% dod = hmrR_ShortChannelRegression_HC(dod, probe, mlActAuto, tIncAuto, rhoSD_ssThresh, PCA, method)
%
% UI NAME:
% shortchannel_regression
%
% DESCRIPTION:
% Regresses out the short channels (with a separation<ssThresh) with an OLS
% or iWLS method. You can chose to only use the first 4 components of PCA analysis
% or just use all the channels.
% This function was created based on ft_nirs_referencechannelsubstraction
% (Fieldtrip toolbox, Artinis external toolbox) and adapted by Helena
% Cockx.
%
% INPUTS:
% dod - SNIRF object containing time course dod (nTpts x nChannels )
% probe - SNIRF object describing the probe - optode positions and wavelengths.
% mlActAuto - list of active channels of the MeasList
% tIncAuto - list of the active time points
% rhoSD_ssThresh - max distance for a short separation measurement. Set =0
%          if you do not want to regress the short separation measurements.
%          Follows the static estimate procedure described in Gagnon et al (2011).
%          NeuroImage, 56(3), 1362?1371.
% PCA - true or false to only use the first 4 principle componenents of the
% short channels or use them all.
% method - 1 (iWLS) or 2 (OLS) to solve the GLM with an iteratively weighted least
% squares model or a ordinary least square method
%
% OUTPUTS:
% dod - SNIRF data containing the time course dod after short channel
% regression
%
% USAGE OPTIONS:
% shortchannel_regression: dod = hmrR_ShortChannelRegression_HC(dod, probe, mlActAuto, tIncAuto, rhoSD_ssThresh, PCA, method)
%
% PARAMETERS:
% rhoSD_ssThresh: 15.0
% PCA: 0
% global_local: 0
% method: 1
%
% TO DO:
% - include the tIncAuto vector
% - make it also work if the input data exist of multiple blocks (for
% iBlk=1:...)
%
function dod = hmrR_ShortChannelRegression_HC(dod, probe, mlActAuto, tIncAuto, rhoSD_ssThresh, PCA, method)
ML = dod.GetMeasList();
ml = dod.GetMeasListSrcDetPairs();
SrcPos = probe.GetSrcPos();
DetPos = probe.GetDetPos();
mlAct = mlActAuto{1};

% find the short channels
lst = 1:size(ml,1);
rhoSD = zeros(length(lst),1);
posM = zeros(length(lst),3);
for iML = 1:length(lst)
  rhoSD(iML) = sum((SrcPos(ml(lst(iML),1),:) - DetPos(ml(lst(iML),2),:)).^2).^0.5;
  posM(iML,:) = (SrcPos(ml(lst(iML),1),:) + DetPos(ml(lst(iML),2),:)) / 2;
end
lstSS = lst(find(rhoSD<=rhoSD_ssThresh & mlAct(lst)==1));
[~,idxSS] = find(ML(:,1)'==ml(lstSS,1) & ML(:,2)'==ml(lstSS, 2));

% mean detrend of the data
dod_temp = bsxfun(@minus, dod.dataTimeSeries, mean(dod.dataTimeSeries,1));
if PCA
  [coeff,score,latent,tsquared,explained] = pca(dod_temp(:, idxSS));
  fprintf('The first 4 PCA components of the short channels contained %.02f%% of the variance. \n',sum(explained(1:4)))
  x = score(:,1:4); % if first 4 components
else
  x = dod_temp(:, idxSS);
end

signal = nan(size(dod_temp));
% reference channel substraction
for i=1:length(dod.measurementList)
  if mlAct(i)~=1
    signal(:,i) = dod_temp(:,i);
    continue
  end
  y = dod_temp(:,i);
  
  if method ==1
    [dmoco, beta, tstat, pval, sigma, CovB, dfe, w, P, f] = ar_glm_final(y, x);
    yhat = x*beta;
    res = y - yhat;
  elseif method==2
    x2 = [repmat(1, size(x, 1), 1) x]; % + constant factor
    beta = x2\y;
    yhat = x2*beta;
    res  = y - yhat;
    beta(1) = [];
  end
  
  signal(:,i) = res;
  
%   % sanity check
%   fprintf('Found the following meaningful (beta>0.5) shallow channels for channel S%d-D%d (wl idx %d):', ML(i,1), ML(i,2), ML(i,4));
%   shIdx = find(beta>0.5);
%   for s=1:numel(shIdx)
%     fprintf('\n\tS%d-D%d (wl idx %d)', ML(idxSS(shIdx(s)),1), ML(idxSS(shIdx(s)),2), ML(idxSS(shIdx(s)),4))
%   end
%   fprintf('\n\n')
  
  % update data
  dod.dataTimeSeries = signal;
  
end
