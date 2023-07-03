close all, clear all
disp('Running...')

%% Define units
global kg meter sec N
kg = 1;
meter = 1;
sec = 1;
N = 1 *(kg*meter/sec^2);

%% Define physical parameters
w               = 0.5 *(meter);                                               % the width of cell
L               = 1.5 *(meter);                                               % the length of cell
rho_lattice     = 30 *(meter^-1);                                           % the density of cross-link
rho_loop        = 20 *(meter^-1);                                           % the density of primary link
l0              = 0.02 *(meter);                                            % the natural length
stiffness       = 5 *(N/sec);                                               % the stiffness of cross-link 
stiffness2      = 10 *(N/sec);                                              % the stiffness of primary link               
damping         = 10 *(N*sec/meter);                                          
mass            = 1 *(kg);                                 
dt              = 0.001 *(sec);                                             % the simulation time step

%% Build nodes
nodes = buildNodes(w, L, rho_lattice, rho_loop);

%% Define GUI
[rows,cols] = size(nodes.mask_left);

S.f = figure('WindowState','maximized');
S.a = axes;
S.h = plot(zeros(2,sum(sum(~isnan(nodes.mask_up)))),zeros(2,sum(sum(~isnan(nodes.mask_up)))),'bo-',...
           zeros(2,sum(sum(isnan(nodes.mask_left)))),zeros(2,sum(sum(isnan(nodes.mask_left)))),'r-');
S.mText = uicontrol('style','text');
xlim([0,L-w]), ylim([0,w*pi])
xlabel('X [meter]', 'fontweight', 'bold');
ylabel('Y [meter]', 'fontweight', 'bold');

drawNodes(S ,nodes, 0);

%% Main
% recording setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v = VideoWriter('relaxation.avi');
% v.FrameRate = 33;
% open(v)
% 
% frame = getframe(gcf);
% writeVideo(v, frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulation part
tic
for i = 1 : 200000
    nodes = updateNode(nodes, mass, l0, stiffness, stiffness2, damping, w, dt);
    
    % monitoring
    if (mod(i,500)==0)
        drawNodes(S ,nodes, i*dt);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         frame = getframe(gcf);
%         writeVideo(v, frame);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close(v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
disp('Done.')

