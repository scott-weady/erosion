%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for generating Fig. 6c,d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
fig = journal_figure([3.25 3],2); %setup figure
load('data/6cd') %load data

dt = 0.001; %time step from simulation
tf = 172; %final time for sigma = 1 case

iexposed = 30/dt; %index of inclusion exposed
ipinch = 80/dt; %index of pinch
texposed = Vn100(iexposed,1)*[1 1]; %for plotting
yexposed = [0 Vn100(iexposed,2)]; % " "
tpinch = Vn100(ipinch,1)*[1 1]; % " "
ypinch = [0 Vn100(ipinch,2)]; % " "

% Define colors
c1 = [0.9363 0.5851 0.3726]; %sigma = 1
c100 = [0.6225 0.3891 0.2478]; %sigma = 100

% Normalize by initial volume
Vol1(:,2) = Vol1(:,2)/Vol1(1,2); 
Vol100(:,2) = Vol100(:,2)/Vol100(1,2); 

% Scale by final time
Vol1(:,1) = Vol1(:,1)/tf; Vol100(:,1) = Vol100(:,1)/tf;
Vn1(:,1) = Vn1(:,1)/tf; Vn100(:,1) = Vn100(:,1)/tf;
texposed = texposed/tf;
tpinch = tpinch/tf;

% Take subset of points for plotting
sp = 1000; %spacing between points for downsampled plotting
Vol1 = Vol1(1:sp:end,:); Vol100 = Vol100(1:sp:end,:);
Vn1 = Vn1(1:sp:end,:); Vn100 = Vn100(1:sp:end,:);
iexposed = iexposed/sp+1;
ipinch = ipinch/sp+1;

%%%% Render volume %%%%
sp1 = subplot(2,1,1);
t = linspace(0,1);
fill([texposed(1) tpinch(1) tpinch(1) texposed(1) texposed(1)],[0 0 1 1 0],c100,'FaceAlpha',0.05,'EdgeColor','none','HandleVisibility','off'), hold on
plot(Vol1(:,1),Vol1(:,2),'Color',c1,'DisplayName','$\sigma = 1$')
plot(Vol100(:,1),Vol100(:,2),'Color',c100,'DisplayName','$\sigma = 100$')
plot(t,(1-t/t(end)).^2,'k--','DisplayName','$(1 - t/t_f)^2$')
legend('location','northeast','edgecolor','none','fontsize',16)
ylabel('volume, $V/V_0$','FontSize',16)
xticklabels({})
xlim([0 1])

plot(texposed,[0 1],'--','Color',0.2*[1 1 1],'LineWidth',1,'HandleVisibility','off')
plot(tpinch,[0 1],'--','Color',0.2*[1 1 1],'LineWidth',1,'HandleVisibility','off')

%%%% Normal velocity %%%%
sp2 = subplot(2,1,2);
fill([texposed(1) tpinch(1) tpinch(1) texposed(1) texposed(1)],[0 0 1 1 0],c100,'FaceAlpha',0.05,'EdgeColor','none'), hold on
plot(Vn1(:,1),Vn1(:,2),'Color',c1)
plot(Vn100(1:ipinch,1),Vn100(1:ipinch,2),'Color',c100)
plot(Vn100(ipinch+1:end,1),Vn100(ipinch+1:end,2),':','Color',c100)
plot(texposed,[0 1],'--','Color',0.2*[1 1 1],'LineWidth',1)
plot(tpinch,[0 1],'--','Color',0.2*[1 1 1],'LineWidth',1)
xlim([0 1]), ylim([0 0.4])
xlabel('time, $t/t_f$','FontSize',16,'interpreter','latex');
ylabel('maximum erosion rate','FontSize',16);

%%%% Format %%%
sp1.Units = 'inches';
sp2.Units = 'inches';

sp2.Position(2) = 0.775;
sp2.Position(3) = 0.85*fig.PaperSize(1);
sp2.Position(4) = sp2.Position(3)*0.4375;
sp2.Position(1) = 0.25+(fig.PaperSize(1)-sp2.Position(3))/2;
sp1.Position(2) = sum(sp2.Position([2 4]))+0.2;
sp1.Position([1 3 4]) = sp2.Position([1 3 4]);
