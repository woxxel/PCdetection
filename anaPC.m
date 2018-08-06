%%%   written by A.Schmidt, last reviewed on August, 6th, 2018

function [PC_fields] = anaPC(activity,bh,para)
  
  if nargin < 3
    para = struct;
    para.f = 15;
    para.nbin = 80;
    para.repnum = 1000;
    para.sigma = 2;
    para.offset = 10;
  end
  PC_fields = struct;
  PC_fields.fields = struct;
  PC_fields.MI = struct;
  
  %%% shuffle only within longrunperiods
  binpos = bh.binpos(bh.longrunperiod);
  activity_lr = activity(bh.longrunperiod);                                            %% only activity from actual times
  spike_times = find(activity_lr);
  spikes = activity_lr(spike_times);
  ISI = diff(spike_times);
  T = length(activity_lr);
  
  PC_fields.firingmap = get_firingmap(binpos(spike_times),spikes,bh.dwelltime);
  [PC_fields.MI.value, PC_fields.MI.binMI] = get_MI(PC_fields.firingmap,bh.norm_dwelltime,para);
  
  MI_rand=zeros(1,para.repnum);
  for L=1:para.repnum
      shuffled_spike_train = spike_shuffling('dithershift',true,spike_times,spikes,T,ISI,2*para.f);
      spike_times_L = find(shuffled_spike_train);
      spikes_L = shuffled_spike_train(spike_times_L);
      
      firingmap = get_firingmap(binpos(spike_times_L),spikes_L,bh.dwelltime);
      [MI_rand(L), ~] = get_MI(firingmap,bh.norm_dwelltime,para);
  end

  PC_fields.MI.frac = PC_fields.MI.value / prctile(MI_rand,95);
  PC_fields.MI.dist = MI_rand;

  status = PC_fields.MI.frac > 1;
  
  if status
    %% determine areas of coding as regions of size > 4 (5cm) of enhanced MI
    %% only consider a cell to be a place cell, if it has a proper firing field
    
    smoothed_MI = imgaussfilt(PC_fields.MI.binMI,para.sigma,'Padding','circular','FilterDomain','spatial');
    periodic_MI = [smoothed_MI(end-(para.offset-1):end) smoothed_MI smoothed_MI(1:para.offset)];
    
    fields = bwareaopen(periodic_MI > 2*nanmean(smoothed_MI),4);
    [bw_comp,~] = bwlabel(fields);
    if bw_comp(1)
      bw_comp(bw_comp==bw_comp(1)) = 0;
    end
    if bw_comp(end)
      bw_comp(bw_comp==bw_comp(end)) = 0;
    end
    
    center_tmp = [];
    for i = unique(bw_comp)
      if i == 0
        continue
      end
      
      loc = find(bw_comp==i);
      centr = mod(round(dot(loc-para.offset,periodic_MI(loc)/sum(periodic_MI(loc)))),para.nbin);
      if ~ismember(centr,center_tmp)
        center_tmp = [center_tmp centr];
      end
    end
    n = length(center_tmp);
    
    PC_fields.status = status && n;
    
    %%% if PC, store some stuff
    if PC_fields.status
      PC_fields.fields.map = fields(para.offset+1:end-para.offset).*smoothed_MI;
      PC_fields.fields.num = n;
      PC_fields.fields.center = center_tmp;
      
%        status_ct = status_ct + 1;
    else
      PC_fields.fields.num = 0;
    end
    
  else
    PC_fields.status = false;
    PC_fields.fields.num = 0;
  end
  
  [PC_fields.max_fr,PC_fields.max_pos] = max(PC_fields.firingmap);
  
end



function [firingmap] = get_firingmap(bins,spikes,dwelltime)
%%% calculates the firing map
  firingmap = zeros(size(dwelltime));
  for i=1:length(bins)
    firingmap(bins(i)) = firingmap(bins(i))+spikes(i);
  end

  firingmap = firingmap./dwelltime;
  firingmap(~dwelltime) = 0;
  
end



function [MI, MI_arr] = get_MI(firingmap,dwelltime,para)
%%% calculates the mutual information according to Skaggs formula, given the provided firingmap and dwelltime
  
  MI_arr = zeros(1,para.nbin);
  mean_firingmap = mean(firingmap);
  
  for k=1:para.nbin
    if ~dwelltime(k)
      MI_arr(k) = NaN;
    elseif firingmap(k)>0
      MI_arr(k) = dwelltime(k) * (firingmap(k)/mean_firingmap)...
            * log2(firingmap(k)/mean_firingmap);
    end
  end
  MI = nansum(MI_arr);
end