% SYNTAX:
% mlActAuto = hmrR_PruneChannelsSCI_HC(dod, probe, mlActMan, tInc, threshold)
%
% UI NAME:
% Prune_ChannelsSCI
%
% DESCRIPTION:
% Prune channels from the measurement list if they have scalp-coupling index 
% lower than the set threshold. This function
% updates MeasListAct based on whether dod 'd' meets these conditions
% as specified by 'threshold'.
% This function was created based on ft_nirs_scalpcouplingindex
% (Fieldtrip toolbox, Artinis external toolbox) and adapted by Helena
% Cockx.
%
% INPUTS:
% d - SNIRF object containing time course dod (nTpts x nChannels )
% probe - SNIRF object describing the probe - optode positions and wavelengths.
% mlActMan - list of active channels of the MeasList
% tInc - list of the active time points 
% threshold - SCI threshold used to prune the channels (default = 0.75)
%
% OUTPUTS:
% mlAct - cell array of all dod blocks - each dod block is an array
%         of true/false for all channels in the contanining dod block
%         specifying active/inactive status. (# of dod blocks x # of Channels)
%
% USAGE OPTIONS:
% Prune_ChannelsSCI: mlActAuto = hmrR_PruneChannelsSCI_HC(dod, probe, mlActMan, tIncMan, threshold)
%
% PARAMETERS:
% threshold: 0.75
%
% TO DO:
% resample data first
%
function mlActAuto = hmrR_PruneChannelsSCI_HC(dod, probe, mlActMan, tIncMan, threshold)

% Init output 
mlActAuto = cell(length(dod),1);

% Check input args
if nargin<5
    disp( 'USAGE: hmrR_PruneChannels(dod, probe, mlActMan, tIncMan, threshold)' )
    return
end
if isempty(tIncMan)
    tIncMan = cell(length(dod),1);
end
if isempty(mlActMan)
    mlActMan = cell(length(dod),1);
end

for iBlk=1:length(dod)
    
    d        = dod(iBlk).GetDataTimeSeries();
    t        = dod(iBlk).GetTime();
    MeasList = dod(iBlk).GetMeasList();
    Lambda   = probe.GetWls();
%     SrcPos   = probe.GetSrcPos();
%     DetPos   = probe.GetDetPos();    
    if isempty(mlActMan{iBlk})
        mlActMan{iBlk} = ones(size(MeasList,1),1);
    end    
    MeasListAct = mlActMan{iBlk};
    MeasListActList = MeasList(find(MeasListAct),:);
    if isempty(tIncMan{iBlk})
        tIncMan{iBlk} = ones(length(t),1);
    end
    tInc = tIncMan{iBlk};
        
    lstInc = find(tInc==1);
    d = d(lstInc,:);
    
    % convert t to fs (assume t is a time vector if length>1)
    if length(t)>1
        fs = 1/(t(2)-t(1));
    end
%     
%     % if sample rate is too high (>10Hz) or we hit Nyquist then resample
%     dataf = d;
%     if fs > 10 || fs <= 5
%       rescfg = [];
%       rescfg.resamplefs = 10;
%       dataf = ft_resampledata(rescfg, dataf);
%     end
    
    % filter the signal in the range of interest
    lpf = 2.4;
    hpf = 0.5;
    % low pass filter
    lpf_norm = lpf / (fs / 2);
    FilterOrder = 3;
    [z, p, k] = butter(FilterOrder, lpf_norm, 'low');
    [sos, g] = zp2sos(z, p, k);
    d2 = filtfilt(sos, g, double(d));
    % high pass filter
    hpf_norm = hpf / (fs / 2);
    FilterOrder = 5;
    [z, p, k] = butter(FilterOrder, hpf_norm, 'high');
    [sos, g] = zp2sos(z, p, k);
    dataf = filtfilt(sos, g, d2);
    
    % loop through channels (assume that related channels are following
    % each other and calculate correlation)
    nChans=size(MeasList,1);
    i=1;
    skipChan = ~MeasListAct;
    sci = ones(nChans, 1);
    while i<nChans
      if skipChan(i)
        i = i+1;
      else
        % find all matching channels
        chanidx=find(MeasList(:,1)==MeasList(i,1) & MeasList(:,2)==MeasList(i,2));
        for ca=1:numel(chanidx)-1
          cidxA = chanidx(ca);
          if skipChan(cidxA)
            continue;
          end
          sci(cidxA) = 0; % initialize this channel
          for cb = ca+1:numel(chanidx) % next channel pair
            cidxB = chanidx(cb);
            if skipChan(cidxB)
              continue;
            end
            sci(cidxB) = 0; % initialize this channel
            
            % correlate the data
            r		= corrcoef(dataf(:,cidxA), dataf(:,cidxB));
            rsum		= r(2);
            
            % sum up the sci, divided by number of 'companion wavelengths'
            sci(cidxA) = sci(cidxA) + rsum/(numel(chanidx)-1);
            sci(cidxB) = sci(cidxB) + rsum/(numel(chanidx)-1);
          end % end for:cb
          
          % mark this channel as 'done'
          skipChan(cidxA) = true;
        end % end for:ca
        
        % skip the last channel as it's either 'done' or has no 'companion'
        skipChan(chanidx(end)) = true;
        
      end % end if:skipChan
    end
    
      % mark the bad channels on the original datain 
  newMeasListAct = (sci>threshold & MeasListAct);  
  fprintf('Marked %d channels as bad. \n', sum(MeasListAct)-sum(newMeasListAct))
  mlActAuto{iBlk} = newMeasListAct;
end
    

