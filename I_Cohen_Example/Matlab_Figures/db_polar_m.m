function [hpol,radlab_h] = db_polar_m(theta, rho_l, rlabshift, clear, linew)
% db_polar_m: to plot the beampattern in dB scale with polar coordinates
%
% Usage: 
%        [hpol,radlab_h] = db_polar_m(theta, rho_l, rlabshift, clear) 
%
%    Input variables:
%       theta:      angles from 0 to 360 degree 
%       rho_l:      beampattern in linear scale
%       rlabshift:  shift the radial lables for impression of negative dB
%                   valudes
%       clear:      make polar background clear: 1 for clear and 0 for
%                   white
%
%
%    Outputs:
%          
% ------
%

if nargin < 1
    error('Requires at least 2 input arguments.');
end;

if ~isequal(size(theta, 1),size(rho_l, 1))
    error('Theta and Rho_l must be the same size.');
end

% Dealing with additional arguments
if nargin < 3,  
    rlabshift = 0;    
end
if nargin < 4,  
    clear = 0;    
end

if isstr(theta) | isstr(rho_l) | isstr(rlabshift)
    error('Input arguments must be numeric.');
end

line_style = 'auto';
labelrotate = 0;

%express the gain into dB stytle between -50 and 0 dB
rho = zeros(size(rho_l));
for i = 1 : size(rho_l, 2),
    rho(:,i) = 20*log10(rho_l(:,i));
    rho(find(rho(:,i) < rlabshift),i) = rlabshift;
    %rho(:,i) =(rho(:,i) + abs(min(rho(:,i))))/abs(min(rho(:,i)));
    rho(:,i) =(rho(:,i) + abs(rlabshift))/abs(rlabshift);
end;

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold on;
    maxrho = max(abs(rho(:)));
    hhh=plot([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho]);
    set(gca,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
    if rticks < 3,
        rticks = 5;
    end;
    if rticks > 5   % see if we can reduce the number
        if rem(rticks,2) == 0
            rticks = rticks/2;
        elseif rem(rticks,3) == 0
            rticks = rticks/3;
        end
    end

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if (~isstr(get(cax,'color')) & ~clear),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end

% draw radial circles
    c82 = cos(82*pi/180);
    s82 = sin(82*pi/180);
    rinc = (rmax-rmin)/rticks;
    cnt = 1;
    for i=(rmin+rinc):rinc:rmax
        hhh = plot(xunit*i,yunit*i,ls,'color',tc,'linewidth',linew,...
                   'handlevisibility','off');
        radlab_h(cnt) = text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
            ['  ' num2str(i*50+rlabshift) ' dB'],'verticalalignment','bottom',...
            'handlevisibility','off');
        cnt = cnt+1;
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = (1:6)*2*pi/12;
    cst = cos(th); snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    plot(rmax*cs,rmax*sn,ls,'color',tc,'linewidth',linew,...
         'handlevisibility','off')

% annotate spokes in degrees
    rt = 1.1*rmax;
    for i = 1:length(th)
        angTextNum = (i*30)+labelrotate    ;
        text(rt*cst(i),rt*snt(i),int2str(angTextNum),...
             'horizontalalignment','center',...
             'handlevisibility','off');
        if i == length(th)
            loc = int2str(0+labelrotate);
        else
            angTextNum = 180+(i*30)+labelrotate;
            %if angTextNum > 180, angTextNum = angTextNum - 360; end
            if angTextNum > 180, 
                angTextNum = angTextNum; 
            end
            loc = int2str(angTextNum);
        end
        text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
             'handlevisibility','off')
    end

% set view to 2-D
    view(2);
% set axis limits
    axis(rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

% transform data to Cartesian coordinates.
colorc = ['b- '; 'r-.'; 'm: '; 'k--'; 'g- '; 'c: '; 'y: '];
if (size(rho, 2) > length(colorc))
    error('Too many different spacing to plot!');
end;
% plot all the beampatterns
for i =1: size(rho, 2)
   xx = rho(:,i).*cos(theta);
   yy = rho(:,i).*sin(theta);   
   % plot data on top of grid
   q = plot(xx, yy, colorc(i, :), 'LineWidth', linew);
end;

if nargout > 0
    hpol = q;
end
if ~hold_state
    set(gca,'dataaspectratio',[1 1 1]), axis off; set(cax,'NextPlot',next);
end
set(get(gca,'xlabel'),'visible','on')
set(get(gca,'ylabel'),'visible','on')
