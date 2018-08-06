%%%   written by A.Schmidt, last reviewed on August, 6th, 2018

function PC_fields = anaPC_frame(s)

pathMouse = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/884';
pathSession = pathcat(pathMouse,sprintf('Session%02d',s));
pathBH = dir(pathcat(pathSession,'*aligned.mat'));
pathBH = pathcat(pathSession,pathBH.name);
bh = load(pathBH);
bh = bh.alignedData.resampled;

pathAct = pathcat(pathSession,'resultsCNMF_MF1_LK1.mat');
act = load(pathAct,'C2','S2');
S = act.S2;
C = act.C2;

nCells = size(C,1);

para = struct;
para.prc = 20;
para.nsd = 3;
para.nbin = 80;
para.binwidth = 1600/para.nbin;

para.f = 15;
para.repnum = 1000;
para.sigma = 2;
para.offset = 10;

PC_fields = struct('fields',struct,'MI',struct,'firingmap',[],'status',cell(nCells,1),'max_fr',cell(nCells,1),'max_pos',cell(nCells,1));
PC_fields
plt = false;
status = 0;
for n = 1:nCells
  
  if ~mod(n,100)
    disp(sprintf('n = %d',n))
  end
  close all
  modeS = prctile(S(n,S(n,:)>0),para.prc);                           %% get mode from overall activity
  activity = floor(sqrt(S(n,:)/(modeS*para.nsd)));
  
  PC_fields(n) = anaPC(activity,bh,para);
  
  if PC_fields(n).status
    status = status + 1;
%      plot_activity(PC_fields(n),C(n,:),S(n,:),bh,para,false)
    
%      waitforbuttonpress
  end
end

disp('Neurons showing receptive fields:')
status



