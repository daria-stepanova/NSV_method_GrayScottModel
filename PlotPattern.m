clear;

Omega = 25000;
font_size=15;
save_fig = 0;

data_u = importdata('Outputs/output_u.txt');
data_v = importdata('Outputs/output_v.txt');

N = size(data_u,2)-1;
h = 1.0/N;

time = data_u(:,end);
lattice = linspace(0,1,N);

f1=figure(1);
clf;
fig = gcf;
fig.Units = 'inches';
fig.Position = [0 3 13.4 7];

ax(1)=subplot(1,2,1);
imagesc(time,lattice,data_u(:,1:N)'./(Omega*h))
set(gca,'YDir','normal')
cu=colorbar;
cu.Label.String = 'u';
cu.Label.Rotation = 0;
cu.Label.FontSize = 20;
cu.Label.Position = [0.4875    1.0743         0];
ylabel('x coordinate','FontSize',font_size)
xlabel('time','FontSize',font_size)
set(gca,'fontsize',font_size,'Xcolor','black','Ycolor','black')

ax(2)=subplot(1,2,2);
imagesc(time,lattice,data_v(:,1:N)'./(Omega*h))
set(gca,'YDir','normal')
cv=colorbar;
cv.Label.String = 'v';
cv.Label.Rotation = 0;
cv.Label.FontSize = 20;
cv.Label.Position = [ 0.4250    1.8596         0];
ylabel('x coordinate','FontSize',font_size)
xlabel('time','FontSize',font_size)
set(gca,'fontsize',font_size,'Xcolor','black','Ycolor','black')
set(gcf,'color','w')

if save_fig==1
    export_fig('TemporalEvolution.png','-transparent','-m6');
end

f2=figure(2);
clf;
fig = gcf;
fig.Units = 'inches';
fig.Position = [0 6 15 3.5];

stairs(lattice,data_u(end,1:N)./(Omega*h),'LineWidth',3,'Color',[0 0.7 0])
hold on
stairs(lattice,data_v(end,1:N)./(Omega*h),'LineWidth',3,'Color',[0.7 0 0])
ylim([0 1.6])
xlabel('x coordinate','FontSize',font_size)
ylabel('concentration','FontSize',font_size)
set(gca,'fontsize',font_size,'Xcolor','black','Ycolor','black')
legend('u','v','FontSize',font_size,'Location','north')
legend boxoff
set(gcf,'color','w')

if save_fig==1
    export_fig('FinalPattern.png','-transparent','-m10');
end

% to plot all the frames
% for ti = 1:size(data_u,1)
%     f2=figure(2);
%     clf;
%     fig = gcf;
%     fig.Units = 'inches';
%     fig.Position = [0 6 15 3.5];
%
%     stairs(lattice,data_u(ti,1:N)./(Omega*h),'LineWidth',3,'Color',[0 0.7 0])
%     hold on
%     stairs(lattice,data_v(ti,1:N)./(Omega*h),'LineWidth',3,'Color',[0.7 0 0])
%     ylim([0 1.6])
%     xlabel('x coordinate','FontSize',font_size)
%     ylabel('concentration','FontSize',font_size)
%     set(gca,'fontsize',font_size,'Xcolor','black','Ycolor','black')
%     legend('u','v','FontSize',font_size,'Location','north')
%     legend boxoff
%     set(gcf,'color','w')
%     export_fig(sprintf('TemporalEvolution_%03d.png',ti),'-m4');
%
% end

