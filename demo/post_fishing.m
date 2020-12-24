
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

%% unconstrained results
fprintf('==============================================================\n')
filename = 'fishing_unc_traj';
fprintf(filename);
data = csvread( [filename,'.csv'] );
fprintf('\n');

t = data(:,1);
xu1 = data(:,2:3);
uu1 = data(:,4);
xu2 = data(:,5:6);
uu2 = data(:,7);
xu3 = data(:,8:9);
uu3 = data(:,10);

figure
subplot(3,1,1)
hold on, grid on, box on
plot(t,xu1(:,1),'-','LineWidth',lnwd,'Color',colors(1,:))
plot(t,xu2(:,1),'--','LineWidth',lnwd,'Color',colors(2,:))
plot(t,xu3(:,1),':','LineWidth',lnwd,'Color',colors(3,:))
ylabel('State x_1')
subplot(3,1,2)
hold on, grid on, box on
plot(t,xu1(:,2),'-','LineWidth',lnwd,'Color',colors(1,:))
plot(t,xu2(:,2),'--','LineWidth',lnwd,'Color',colors(2,:))
plot(t,xu3(:,2),':','LineWidth',lnwd,'Color',colors(3,:))
ylabel('State x_2')
subplot(3,1,3)
hold on, grid on, box on
plot(t,uu1,'-','LineWidth',lnwd,'DisplayName','sigma = 0.0','Color',colors(1,:))
plot(t,uu2,'--','LineWidth',lnwd,'DisplayName','sigma = 0.1','Color',colors(2,:))
plot(t,uu3,':','LineWidth',lnwd,'DisplayName','sigma = 1.0','Color',colors(3,:))
legend('show')
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

%% constrained results
fprintf('==============================================================\n')
filename = 'fishing_con_traj';
fprintf(filename);
data = csvread( [filename,'.csv'] );
fprintf('\n');

t = data(:,1);
xunc = data(:,2:3);
uunc = data(:,4);
xcon = data(:,5:6);
ucon = data(:,7);

figure
subplot(3,1,1)
hold on, grid on, box on
plot(t,xunc(:,1),'-','LineWidth',lnwd,'Color',colors(1,:))
plot(t,xcon(:,1),'--','LineWidth',lnwd,'Color',colors(2,:))
ylabel('State x_1')
subplot(3,1,2)
hold on, grid on, box on
plot(t,xunc(:,2),'-','LineWidth',lnwd,'Color',colors(1,:))
plot(t,xcon(:,2),'--','LineWidth',lnwd,'Color',colors(2,:))
ylabel('State x_2')
subplot(3,1,3)
hold on, grid on, box on
plot(t,uunc,'-','LineWidth',lnwd,'DisplayName','unconstrained','Color',colors(1,:))
plot(t,ucon,'--','LineWidth',lnwd,'DisplayName','constrained','Color',colors(2,:))
legend('show')
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



