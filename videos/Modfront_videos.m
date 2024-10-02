%% This script generates the videos of different pattern interfaces
% Note: If you are using a Matlab-version prior to 2022a, you have to
% change "clim" to "caxis" below.


%% set up 2d mesh
Nx = 500;
Ny = 100;
xmin = -10;
xmax = 10;
ymax = 2;
ymin = -2;
x = linspace(xmin,xmax,Nx);
y = linspace(ymin,ymax,Ny);

[X,Y] = meshgrid(x,y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H_{1,+} -> T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = @(t,x,y) (1/2)*(tanh(-0.5*(x-t))+1).*(cos(5*x)+cos(5/2*(-x+sqrt(3)*y))+cos(-5/2*(x+sqrt(3)*y)))-2;
f = @(t,x,y) (1./sqrt(3.*exp(x-t)+1)).*(cos(5*x)+cos(5/2*(-x+sqrt(3)*y))+cos(-5/2*(x+sqrt(3)*y)));

% maximal amplitude
maxAmp = max(abs(f(100,X,Y)),[],"all");

heights = linspace(min(f(0,X,Y),[],"all"),max(f(0,X,Y),[],"all"),9);
% contourf(X,Y,f(0,X,Y),heights,'k')
% pbaspect([xmax/ymax 1 1])

v = VideoWriter("Modfronts-H-to-T.mp4","MPEG-4");
v.FrameRate = 30;
v.Quality = 100;
open(v)

fig = figure('units','pixels','position',[100 200 2000 400]);
axes(fig,'units','pixels','position',[0 0 2000 400]);

t = linspace(-15,15,400);

for ii = 1:numel(t)
    % hold on
    s = surfc(X,Y,f(t(ii),X,Y));
    axis off;
    view(2);
    pbaspect([xmax/ymax 1 1])
    grid off;
    s(1).EdgeColor="none";
    colormap jet
    % clim([min(f(0,X,Y),[],"all"),max(f(0,X,Y),[],"all")]);
    clim([-maxAmp,maxAmp]);
    ax = gca;
    ax.ZLim(2) = max(f(0,X,Y),[],"all");
    s(2).ZLocation = 'zmax';
    s(2).LevelList = heights;
    s(2).EdgeColor="k";
    % contour(X,Y,f(t(ii),X,Y),heights,'k');
    % hold off
    drawnow;

    % contourf(X,Y,f(t(ii),X,Y),heights,'k');
    % pbaspect([xmax/ymax 1 1])
    % clim([min(f(0,X,Y),[],"all"),max(f(0,X,Y),[],"all")]);


    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% R -> T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(t,x,y) (1./sqrt(3.*exp((x-t))+1)).*(cos(5*x));

heights = linspace(min(f(0,X,Y),[],"all"),max(f(0,X,Y),[],"all"),6);

% maximal amplitude
maxAmp = max(abs(f(100,X,Y)),[],"all");

v = VideoWriter("Modfronts-R-to-T.mp4","MPEG-4");
v.FrameRate = 30;
v.Quality = 100;
open(v)

fig = figure('units','pixels','position',[100 200 2000 400]);
axes(fig,'units','pixels','position',[0 0 2000 400]);

t = linspace(-15,15,400);
% t = 0;

for ii = 1:numel(t)
    s = surfc(X,Y,f(t(ii),X,Y));
    axis off;
    view(2);
    pbaspect([xmax/ymax 1 1])
    grid off;
    s(1).EdgeColor="none";
    colormap jet
    % clim([min(f(0,X,Y),[],"all"),max(f(0,X,Y),[],"all")]);
    clim([-maxAmp,maxAmp]);
    ax = gca;
    ax.ZLim(2) = max(f(0,X,Y),[],"all");
    s(2).ZLocation = 'zmax';
    s(2).LevelList = heights(2:end-1);
    s(2).EdgeColor="k";
    drawnow;

    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H_{2,-} -> T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(t,x,y) (1./sqrt(3.*exp(x-t)+1)).*(-cos(5*x)+cos(5/2*(-x+sqrt(3)*y))+cos(-5/2*(x+sqrt(3)*y)));

heights = linspace(min(f(0,X,Y),[],"all"),max(f(0,X,Y),[],"all"),9);

% maximal amplitude
maxAmp = max(abs(f(100,X,Y)),[],"all");

v = VideoWriter("Modfronts-down-H-to-T.mp4","MPEG-4");
v.FrameRate = 30;
v.Quality = 100;
open(v)

fig = figure('units','pixels','position',[0 100 2000 400]);
axes(fig,'units','pixels','position',[0 0 2000 400]);

t = linspace(-15,15,400);

for ii = 1:numel(t)
    s = surfc(X,Y,f(t(ii),X,Y));
    axis off;
    view(2);
    pbaspect([xmax/ymax 1 1])
    grid off;
    s(1).EdgeColor="none";
    colormap jet
    % clim([min(f(0,X,Y),[],"all"),max(f(0,X,Y),[],"all")]);
    clim([-maxAmp,maxAmp]);
    ax = gca;
    ax.ZLim(2) = max(f(0,X,Y),[],"all");
    s(2).ZLocation = 'zmax';
    s(2).LevelList = heights;
    s(2).EdgeColor="k";
    drawnow;

    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H_{1,+} -> H_{2,-}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate heteroclinic orbit
M0 = 1;
N0 = 1;
K0 = -1;
K2 = -2;

% vector field
vecfieldRev = @(t,A) -[-M0*A(1)-N0*A(2).^2-K0*A(1).^3-2*K2*A(1).*A(2).^2;-M0*A(2)-N0*A(1).*A(2)-(K0+K2)*A(2).^3-K2*A(2).*A(1).^2];

roll = [sqrt(-M0/K0),0];
hex = [(-N0-sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2)),(-N0-sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2))];
hexdown = [(-N0+sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2)),-(-N0+sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2))];

eps = 1e-15;

init = hexdown + eps*[2,1];

% calculate heteroclinic H_{2,-} --> H_{1,+} in time-reversed system
[thetRev,xhetRev] = ode45(vecfieldRev,[-50,50],init);

% time-reversal
thet = -flipud(thetRev);
xhet = flipud(xhetRev);

f = @(t,x,y) interp1(thet,xhet(:,1),0.5*(x-t)).*cos(5*x) + interp1(thet,xhet(:,2),0.5*(x-t)).*(cos(5/2*(-x+sqrt(3)*y))+cos(-5/2*(x+sqrt(3)*y)));
maxAmp = max(xhet(:,1)+2*xhet(:,2));

t0 = -20;
heights = linspace(min(f(0,X,Y),[],"all"),max(f(0,X,Y),[],"all"),9);
% contourf(X,Y,f(0,X,Y),heights,'k')
% pbaspect([xmax/ymax 1 1])

v = VideoWriter("Modfronts-H-to-down-H.mp4","MPEG-4");
v.FrameRate = 30;
v.Quality = 100;
open(v)

fig = figure('units','pixels','position',[0 100 2000 400]);
axes(fig,'units','pixels','position',[0 0 2000 400]);

t = linspace(-15,15,400)+t0;

for ii = 1:numel(t)
    s = surfc(X,Y,f(t(ii),X,Y));
    axis off;
    view(2);
    pbaspect([xmax/ymax 1 1])
    grid off;
    s(1).EdgeColor="none";
    colormap jet
    clim([-maxAmp,maxAmp]);
    ax = gca;
    ax.ZLim(2) = max(f(0,X,Y),[],"all");
    s(2).ZLocation = 'zmax';
    s(2).LevelList = heights;
    s(2).EdgeColor="k";
    drawnow;

    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H_{1,+} -> R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate heteroclinic orbit
M0 = 0.5;
N0 = 1;
K0 = -1;
K2 = -2;

% vector field
vecfieldRev = @(t,A) -[-M0*A(1)-N0*A(2).^2-K0*A(1).^3-2*K2*A(1).*A(2).^2;-M0*A(2)-N0*A(1).*A(2)-(K0+K2)*A(2).^3-K2*A(2).*A(1).^2];

roll = [sqrt(-M0/K0),0];
hex = [(-N0-sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2)),(-N0-sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2))];
hexdown = [(-N0+sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2)),-(-N0+sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2))];

eps = 1e-6;

init = roll + eps*[0,1];

% calculate heteroclinic H_{2,-} --> H_{1,+} in time-reversed system
[thetRev,xhetRev] = ode45(vecfieldRev,[-50,50],init);

% time-reversal
thet = -flipud(thetRev);
xhet = flipud(xhetRev);

f = @(t,x,y) interp1(thet,xhet(:,1),2*(x-t),'linear','extrap').*cos(5*x) + interp1(thet,xhet(:,2),2*(x-t),'linear','extrap').*(cos(5/2*(-x+sqrt(3)*y))+cos(-5/2*(x+sqrt(3)*y)));

maxAmp = max(xhet(:,1)+2*xhet(:,2));

t0 = 2.5;
heights = linspace(min(f(0,X,Y),[],"all"),max(f(0,X,Y),[],"all"),9);
% contourf(X,Y,f(0,X,Y),heights,'k')
% pbaspect([xmax/ymax 1 1])

v = VideoWriter("Modfronts-H-to-R.mp4","MPEG-4");
v.FrameRate = 30;
v.Quality = 100;
open(v)

fig = figure('units','pixels','position',[0 100 2000 400]);
axes(fig,'units','pixels','position',[0 0 2000 400]);

t = linspace(-15,15,400)+t0;

for ii = 1:numel(t)
    s = surfc(X,Y,f(t(ii),X,Y));
    axis off;
    view(2);
    pbaspect([xmax/ymax 1 1])
    grid off;
    s(1).EdgeColor="none";
    colormap jet
    clim([-maxAmp,maxAmp]);
    ax = gca;
    ax.ZLim(2) = maxAmp;
    s(2).ZLocation = 'zmax';
    s(2).LevelList = heights;
    s(2).EdgeColor="k";
    drawnow;

    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% two stage invasion: H_{1,+} -> R -> T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate heteroclinic orbit
M0 = 0.5;
N0 = 1;
K0 = -1;
K2 = -2;

% vector field
vecfield = @(t,A) [-M0*A(1)-N0*A(2).^2-K0*A(1).^3-2*K2*A(1).*A(2).^2;-M0*A(2)-N0*A(1).*A(2)-(K0+K2)*A(2).^3-K2*A(2).*A(1).^2];

roll = [sqrt(-M0/K0),0];
hex = [(-N0-sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2)),(-N0-sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2))];
hexdown = [(-N0+sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2)),-(-N0+sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2))];

% select point on orbit
eps = 1e-3;
pt = roll + eps*[-1,1];

% forward orbit
[tfwd,xfwd] = ode45(vecfield,[0,50],pt);

% backward orbit
[tbwd,xbwd] = ode45(@(t,x)-vecfield(t,x),[0,50],pt);

% reverse backward orbit and glue orbits together
tbwdRev = -flipud(tbwd);
xbwdRev = flipud(xbwd);

thet = [tbwdRev(1:end-1);tfwd];
xhet = [xbwdRev(1:end-1,:);xfwd];

% plot(thet,xhet(:,1));

f = @(t,x,y) interp1(thet,xhet(:,1),1*(x-t),'linear','extrap').*cos(5*x) + interp1(thet,xhet(:,2),1*(x-t),'linear','extrap').*(cos(5/2*(-x+sqrt(3)*y))+cos(-5/2*(x+sqrt(3)*y)));

maxAmp = max(xhet(:,1)+2*xhet(:,2));

% heights = linspace(min(f(-50,X,Y),[],"all"),max(f(-50,X,Y),[],"all"),9);
heights = linspace(-maxAmp,maxAmp,8);

v = VideoWriter("Modfronts-H-to-R-to-T.mp4","MPEG-4");
v.FrameRate = 30;
v.Quality = 100;
open(v)

fig = figure('units','pixels','position',[0 100 2000 400]);
axes(fig,'units','pixels','position',[0 0 2000 400]);

t = linspace(-25,45,400);
% t = 45;

for ii = 1:numel(t)
    s = surfc(X,Y,f(t(ii),X,Y));
    axis off;
    view(2);
    pbaspect([xmax/ymax 1 1])
    grid off;
    s(1).EdgeColor="none";
    colormap jet
    clim([-maxAmp,maxAmp]);
    ax = gca;
    ax.ZLim(2) = maxAmp;
    s(2).ZLocation = 'zmax';
    s(2).LevelList = heights;
    s(2).EdgeColor="k";
    drawnow;

    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% two stage invasion: H_{1,+} -> MM -> T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate heteroclinic orbit
M0 = 2;
N0 = 1;
K0 = -1;
K2 = -2;

% vector field
vecfield = @(t,A) [-M0*A(1)-N0*A(2).^2-K0*A(1).^3-2*K2*A(1).*A(2).^2;-M0*A(2)-N0*A(1).*A(2)-(K0+K2)*A(2).^3-K2*A(2).*A(1).^2];

roll = [sqrt(-M0/K0),0];
hex = [(-N0-sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2)),(-N0-sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2))];
hexdown = [(-N0+sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2)),-(-N0+sqrt(N0^2-4*M0*(K0+2*K2)))/(2*(K0+2*K2))];
mm = [N0/(K0-K2),(1/(K0-K2))*sqrt(-(K0*N0^2+(K0-K2)^2*M0)/(K0+K2))];

% select point on orbit
eps = 1e-4;
pt = mm + eps*[-1,0];

% forward orbit
[tfwd,xfwd] = ode45(vecfield,[0,50],pt);

% backward orbit
[tbwd,xbwd] = ode45(@(t,x)-vecfield(t,x),[0,50],pt);

% reverse backward orbit and glue orbits together
tbwdRev = -flipud(tbwd);
xbwdRev = flipud(xbwd);

thet = [tbwdRev(1:end-1);tfwd];
xhet = [xbwdRev(1:end-1,:);xfwd];

% plot(thet,xhet(:,1));

f = @(t,x,y) interp1(thet,xhet(:,1),1*(x-t),'linear','extrap').*cos(5*x) + interp1(thet,xhet(:,2),1*(x-t),'linear','extrap').*(cos(5/2*(-x+sqrt(3)*y))+cos(-5/2*(x+sqrt(3)*y)));

maxAmp = max(xhet(:,1)+2*xhet(:,2));

% heights = linspace(min(f(-50,X,Y),[],"all"),max(f(-50,X,Y),[],"all"),9);
heights = linspace(-maxAmp,maxAmp,8);

v = VideoWriter("Modfronts-H-to-MM-to-T.mp4","MPEG-4");
v.FrameRate = 30;
v.Quality = 100;
open(v)

fig = figure('units','pixels','position',[0 100 2000 400]);
axes(fig,'units','pixels','position',[0 0 2000 400]);

t = linspace(-25,45,400);
% t = 45;

for ii = 1:numel(t)
    s = surfc(X,Y,f(t(ii),X,Y));
    axis off;
    view(2);
    pbaspect([xmax/ymax 1 1])
    grid off;
    s(1).EdgeColor="none";
    colormap jet
    clim([-maxAmp,maxAmp]);
    ax = gca;
    ax.ZLim(2) = maxAmp;
    s(2).ZLocation = 'zmax';
    s(2).LevelList = heights;
    s(2).EdgeColor="k";
    drawnow;

    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)