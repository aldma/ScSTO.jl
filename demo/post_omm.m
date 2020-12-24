
close all
clear all %#ok
clc

addpath('~/Documents/MATLAB/matlab2tikz/src/');
settings.save_figures = true;
lnwd = 1.5;

colors = [252, 141, 98;...
          141, 160, 203;...
          102, 194, 165]/256;
      
if isfield(settings,'save_figures') && settings.save_figures
    fig_style = get_figure_style();
    print_tikz = @(figname) matlab2tikz( [figname,'.tikz'],...
                                      'width', '\columnwidth',...
                                      'height', '\columnwidth',...
                                      'checkForUpdates', false,...
                                      'parseStrings', false,...
                                      'showInfo', false );
    fprintf('prepared for saving figures!\n');
else
    fprintf('NOT saving figures...\n');
end

%% initial guess
fprintf('==============================================================\n')
filename = 'omm_init_traj';
fprintf(filename);
data = csvread( [filename,'.csv'] );
fprintf('\n');

t = data(:,1);
xu1 = data(:,2:3);
uu1 = data(:,4);
xu2 = data(:,5:6);
uu2 = data(:,7);

figure
subplot(3,1,1)
hold on, grid on, box on
plot(t,xu1(:,1),'--','LineWidth',lnwd,'Color',colors(1,:))
plot(t,xu2(:,1),'-','LineWidth',lnwd,'Color',colors(2,:))
ylabel('State x_1')
ylim([0.5,1.1])
subplot(3,1,2)
hold on, grid on, box on
plot(t,xu1(:,2),'--','LineWidth',lnwd,'Color',colors(1,:))
plot(t,xu2(:,2),'-','LineWidth',lnwd,'Color',colors(2,:))
ylabel('State x_2')
subplot(3,1,3)
hold on, grid on, box on
plot(t,uu1,'--','LineWidth',lnwd,'DisplayName','initial guess','Color',colors(1,:))
plot(t,uu2,'-','LineWidth',lnwd,'DisplayName','sigma = 0','Color',colors(2,:))
legend('show')
ylim([-0.1,2.1])
xlabel('Time t')
ylabel('Control w')
drawnow

if isfield(settings,'save_figures') && settings.save_figures
    % apply figure style
    grid on, box on
    hgexport(gcf,'dummy',fig_style,'applystyle', true);
    % store either .tikz  or .fig
    cleanfigure();
    print_tikz( filename );
    pause(0.5)
    close
end

%% switching costs
fprintf('==============================================================\n')
filename = 'omm_swc_traj';
fprintf(filename);
data = csvread( [filename,'.csv'] );
fprintf('\n');

t = data(:,1);
x1 = data(:,2:3);
u1 = data(:,4);
x2 = data(:,5:6);
u2 = data(:,7);
x3 = data(:,8:9);
u3 = data(:,10);

figure
subplot(3,1,1)
hold on, grid on, box on
plot(t,x1(:,1),'-','LineWidth',lnwd,'Color',colors(2,:))
plot(t,x2(:,1),'--','LineWidth',lnwd,'Color',colors(1,:))
plot(t,x3(:,1),':','LineWidth',lnwd,'Color',colors(3,:))
ylabel('State x_1')
ylim([0.5,1.1])
subplot(3,1,2)
hold on, grid on, box on
plot(t,x1(:,2),'-','LineWidth',lnwd,'Color',colors(2,:))
plot(t,x2(:,2),'--','LineWidth',lnwd,'Color',colors(1,:))
plot(t,x3(:,2),':','LineWidth',lnwd,'Color',colors(3,:))
ylabel('State x_2')
subplot(3,1,3)
hold on, grid on, box on
plot(t,u1,'-','LineWidth',lnwd,'Color',colors(2,:),'DisplayName','sigma = 0')
plot(t,u2,'--','LineWidth',lnwd,'Color',colors(1,:),'DisplayName','sigma = 25')
plot(t,u3,':','LineWidth',lnwd,'Color',colors(3,:),'DisplayName','sigma = 30')
legend('show')
ylim([-0.1,2.1])
xlabel('Time t')
ylabel('Control w')
drawnow

if isfield(settings,'save_figures') && settings.save_figures
    % apply figure style
    grid on, box on
    hgexport(gcf,'dummy',fig_style,'applystyle', true);
    cleanfigure();
    print_tikz( filename );
    pause(0.5)
    close
end

% end of file
fprintf('\nThat`s all folks!\n')



