%%% function to shuffle spike trains according to Gansel, 2012
%%%
%%% inputs:
%%%         mode  - specify mode for shuffling
%%%               'shift'       - shift spike train by a fixed offset (default)
%%%                       provide values as (mode,shuffle_peaks,spike_train)
%%%               'dither'      - dither each spike by an independently drawn random value (max = w)
%%%                       provide values as (mode,shuffle_peaks,spike_times,spikes,T,ISI,w,shuffle_spikes)
%%%               'dithershift' - combination of 'shift' and 'dither' method
%%%                       provide values as (mode,shuffle_peaks,spike_times,spikes,T,ISI,w,shuffle_spikes)
%%%               'dsr'         - (dither-shift-reorder), same as 'dithershift' but with random reordering of consecutive ISIs < w (?)
%%%                       provide values as (mode,shuffle_peaks,spike_times,spikes,T,ISI,w,shuffle_spikes)
%%%
%%%         shuffle_peaks  - boolean: should assignment "spikes" to "spike_times" be shuffled?
%%%
%%%         spike_train - spike train as binary array (should be replaced by ISI & T
%%%
%%%         spike_times - frames at which spikes happen
%%%
%%%         spikes      - number of spikes happening at times "spike_times"
%%%
%%%         T     - length of the overall recording (= maximum value for new spike time)
%%%
%%%         ISI   - InterSpike Intervalls of the spike train
%%%
%%%         w     - maximum dithering (~1/(2*rate)?)
%%%
%%% ouputs:
%%%         new_spike_train - shuffled spike train
%%%
%%%   written by A.Schmidt, last reviewed on August, 6th, 2018


function [new_spike_train] = spike_shuffling(mode,shuffle_peaks,varargin)

  switch mode
    case 'shift'
      spike_train = varargin{1};
      [new_spike_train,~] = shift_spikes(spike_train);
      if shuffle_peaks
        spike_times = find(new_spike_train);
        spikes = new_spike_train(spike_times);
        new_spike_train(spike_times) = spikes(randperm(length(spike_times)));        %% shuffle spike numbers
      end
      
    case 'dither'
      
      assert(nargin>=4,'You did not provide enough input. Please check the function description for further information.')
      [spike_times,spikes,T,ISI,w] = get_input_dither(varargin);
      
      new_spike_train = dither_spikes(spike_times,spikes,T,ISI,w,shuffle_peaks);
      
    case 'dithershift'
      
      assert(nargin>=4,'You did not provide enough input. Please check the function description for further information.')
      [spike_times,spikes,T,ISI,w] = get_input_dither(varargin);
      
      new_spike_train = dither_spikes(spike_times,spikes,T,ISI,w,shuffle_peaks);
      [new_spike_train,shift] = shift_spikes(new_spike_train);
      
    case 'dsr'
    
      disp('not yet implemented')
      new_spike_train = NaN;
      
  end
  
  plt = false;
  if plt
    
    if ~exist('spike_train','var')
      spike_train = zeros(1,T);
      spike_train(spike_times) = spikes;
    end
    ISI = get_ISI(spike_train);
    newISI = get_ISI(new_spike_train);
    
    figure('position',[500 500 1200 900])
    subplot(3,1,1)
    plot(spike_train)
    subplot(3,1,2)
    plot(new_spike_train)
    title('new spike train')
    
    subplot(3,2,5)
    hold on
    histogram(log10(ISI),linspace(-2,2,51),'FaceColor','b')
    histogram(log10(newISI),linspace(-2,2,51),'FaceColor','r')
    hold off
    
    waitforbuttonpress;
  end
end



function [new_spike_train,shift] = shift_spikes(spike_train)

  shift = randi(length(spike_train));
  new_spike_train = [spike_train(shift:end),spike_train(1:shift-1)];    %% shift spike train
  
end



function [spike_times,spikes,T,ISI,w] = get_input_dither(input)
  
  if length(input{2}) == 1   
    spike_train = input{1};
    spike_times = find(spike_train);
    spikes = spike_train(spike_times);
    T = length(spike_train);
    ISI = diff(spike_times);
    w = input{2};
  else
    spike_times = input{1};
    spikes = input{2};
    T = input{3};
    ISI = input{4};
    w = input{5};
  end
end


function [new_spike_train] = dither_spikes(spike_times,spikes,T,ISI,w,shuffle_peaks)
  
  nspike_times = length(spike_times);
  
  dither = min(ISI-1,2*w)/2;
  
  r = 2*(rand(1,length(ISI)-1)-0.5);
  
  for i=2:length(ISI)   %% probability of being left or right of initial spike should be equal! (otherwise, it destroys bursts!)
    spike_times(i) = spike_times(i) + min(0,r(i-1))*ISI(i-1) + max(0,r(i-1))*ISI(i);
  end
  spike_times = round(spike_times);
  
  if shuffle_peaks
    spikes = spikes(randperm(nspike_times));
  end
  
  new_spike_train = zeros(1,T);
  for i=1:nspike_times
    t = spike_times(i);
    new_spike_train(t) = new_spike_train(t) + spikes(i);
  end
  
end


function [ISI] = get_ISI(spike_train)
  
  %% this part effectively splits up spike bursts
  spike_times = find(spike_train);
  idx_old = 1;
  new_spike_times = [];
  for t = find(spike_train>1)
    idx_new = find(spike_times==t);
    nspikes = spike_train(t);
    
    new_spike_times = [new_spike_times spike_times(idx_old:idx_new-1) t+linspace(0,1-1/nspikes,nspikes)];
    idx_old = idx_new+1;
  end
  new_spike_times = [new_spike_times spike_times(idx_old:end)];
  ISI = diff(new_spike_times);
  
end