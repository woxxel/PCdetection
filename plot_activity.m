%%%   written by A.Schmidt, last reviewed on August, 6th, 2018

function plot_activity(PC_fields,C,S,bh,para,sv)
  
%    pathSave = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Discussion/';
  
  modeS = prctile(S(S>0),para.prc);                    %% get mode from overall activity
  activity = floor(sqrt(S/(modeS*para.nsd)));         %% only activity from actual times
  
  figure('position',[100 100 2400 800])
  ax1 = axes('position',[0.05,0.725,0.9,0.25]);
  ax2 = axes('position',[0.05,0.475,0.9,0.25]);
  ax3 = axes('position',[0.05,0.1,0.9,0.375]);
  
  ax_MI = axes('position',[0.01,0.5,0.1,0.1]);
  
  hold(ax_MI,'on')
  histogram(ax_MI,PC_fields.MI.dist,'FaceColor','b')
  plot(ax_MI,PC_fields.MI.value,0,'rx','MarkerSize',10)
  text(ax_MI,PC_fields.MI.value-0.2,60,'unshuffled','FontSize',8)
  hold(ax_MI,'off')
  
  set(ax_MI,'YTick',[])
  set(ax_MI,'YColor','none')
  title(ax_MI,'Shuffled MI')
  
  hold(ax2,'on')
  plot(ax2,bh.time,S,'k')
  plot(ax2,[0,bh.time(end)],[modeS,modeS]*para.nsd,'r--','LineWidth',2)
  
  b1 = barh(ax2,0,1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','None','DisplayName','Place field');
  b2 = barh(ax2,0,1,'FaceColor','r','DisplayName','Average firing rate');
  
  hold(ax2,'off')
  ylim(ax2,[0,max(S)*0.8])
  set(ax2,'XTick',[])
  set(ax2,'YAxisLocation','right')
  set(ax2,'ytick',[])
  ylabel(ax2,'deconvolved')
  
  
  
  hold(ax3,'on')
  
  if PC_fields.status > 0
    barh(ax3,(PC_fields.fields.map>0)*bh.time(end),1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','None')
  end
  
  barh(ax3,-PC_fields.firingmap*(bh.time(end)/10)/max(PC_fields.firingmap),'r')
  plot(ax3,bh.time,bh.position/para.binwidth)
  
  if length(activity) > 0
    idx = activity & bh.longrunperiod;
    activity_lr = zeros(1,8989);
    activity_lr(idx) = activity(idx);
    scatter(ax3,bh.time(find(activity_lr>0)),bh.position(activity_lr>0)/para.binwidth,3*activity_lr(activity_lr>0)+5,'r','fill')
    
    idx = activity & ~bh.longrunperiod;
    activity_nlr = zeros(1,8989);
    activity_nlr(idx) = activity(idx);
    scatter(ax3,bh.time(find(activity_nlr>0)),bh.position(activity_nlr>0)/para.binwidth,3*activity_nlr(activity_nlr>0)+5,'k','fill')
  end
  
%    if PC_fields(c).status(s) > 0
%      title(ax3,sprintf('PC, MI = %4.2g (%4.2g)',PC_fields(c).MI(s).value,PC_fields(c).MI(s).frac))
%      text(ax3,0.1,0.5,sprintf('MI = %4.2g',PC_fields(c).MI(s).value))
%    else
%      title(ax3,sprintf('nPC, MI = %4.2g (%4.2g)',PC_fields(c).MI(s).value,PC_fields(c).MI(s).frac))
%    end
  
  hold(ax1,'on')
  plot(ax1,bh.time,C,'k','DisplayName','Calcium data')
  plot(ax1,[0 0],[0 0],'b','DisplayName','Mouse location')
  plot(ax1,[0 0],[0 0],'r--','Linewidth',2,'DisplayName','Activity threshold')
  
  scatter(ax1,bh.time(activity_lr>0),C(activity_lr>0),activity_lr(activity_lr>0)*3+5,'ro','filled','DisplayName','AP during activity')
  scatter(ax1,bh.time(activity_nlr>0),C(activity_nlr>0),activity_nlr(activity_nlr>0)*3+5,'ko','filled','DisplayName','AP during rest')
  hold(ax1,'off')
%    xlim(ax1,[-bh.time(end)/10,bh.time(end)])
  
  set(ax1,'XTick',[])
  set(ax1,'YTick',[])
  set(ax1,'YAxisLocation','right')
%    set(ax1,'ycolor','none')
%    ax1.YAxisLine = 'off';
  ylabel(ax1,'Ca^{2+} activity')
  
  hold(ax3,'off')
  set(ax3,'YAxisLocation','right')
  ylim(ax3,[0 para.nbin])
  
  set(ax3,'YTick',linspace(0,80,5))
  set(ax3,'YTickLabels',linspace(0,para.nbin,5))
  xlabel(ax3,'t [s]')
  ylabel(ax3,'Location [cm]')
%    xlim(ax3,[])
  linkaxes([ax1,ax2,ax3],'x')
  xlim(ax3,[-bh.time(end)/10,bh.time(end)])
  legend(ax1,'Location','NorthWest')
  legend(ax2,[b1,b2],'Location','NorthWest')
%    if sv
%      path = sprintf('%sfromCaToPC.png',pathSave);
%      print(path,'-dpng','-r600')
%      disp(sprintf('saved image in %s',path))
%    end
%  end