%--------------------------------
% Simulation of an eroding valley under a glacier
%   Starts at steady-state glacier shape from original glacier assigment
%   Approximates changes in subglacial water network by changes in 
%       bulk hydraulic conductivity
%   Temperature forcing is handled by oscillating melt input
%--------------------------------
% Franklin Hinckley
% 25 April 2016
%--------------------------------
%
%--------------------------------

%% Clean up workspace
clearvars
close all
clc

%% Constants for glacier [MacGregor et al 2000]
g     = 9.81;    % standard gravity [m/s^2]
A     = 2.1e-16; % Arrhenius Constant [Pa^-3*yr^-1]
gamma = 0.01;    % net mass balance gradient [m/yr/m]
rhoI  = 917;     % density of ice [kg/m^3]
rhoW  = 1000;    % density of water [kg/m^3]
usl0  = 0.0012;  % sliding coefficient [m/yr/Pa]
edot0 = 0.0001;  % erosion coefficient []
bmax  = 3;       % maximum accumulaton [m/yr]

%% Constants for subglacial drainage
k = 4e-5*86400 ;   % initial hydraulic conductivity [m/yr] [Fountain and Walder]
p = 0.05;          % porosity of glacier
scl_flow  = 1e-7;  % scale factor for conductivity as a function of water flow
scl_slide = 1e-9;  % scale factor for conductivity as a function of sliding
scl_close = 1e-14; % scale factor for conductivity as a function of ice thickness

%% Set up initial valley shape
dx = 250;       % [m]
len = 200*1000; % [m]
x = 0:dx:len;

SR = 20/1000;  % slope [m/m]

zMax = 3000; % [m]
zR = zMax - SR*x;

%% Time array
dt   = 1/365/24/4; % [yr]
tSim = 10;         % [yr]
t    = 0:dt:tSim;

%% Initial glacier shape
load Hsteady250
H = Hsteady250'; % steady state glacier from original code

%% Allocate output
% Specify how often to save
saveInd = 1000;

% Allocate output
zR_S  = zeros(length(x)  , floor(length(t)/saveInd));
H_S   = zeros(length(x)  , floor(length(t)/saveInd));
HW_S  = zeros(length(x)  , floor(length(t)/saveInd));
k_S   = zeros(length(x)  , floor(length(t)/saveInd));
usl_S = zeros(length(x)-1, floor(length(t)/saveInd));
Q_S   = zeros(length(x)+1, floor(length(t)/saveInd));
QW_S  = zeros(length(x)  , floor(length(t)/saveInd));
m_S   = zeros(length(x)  , floor(length(t)/saveInd));
ELA_S = zeros(floor(length(t)/saveInd),1);
V     = zeros(floor(length(t)/saveInd),1);
VW    = zeros(floor(length(t)/saveInd),1);
hydro = zeros(floor(length(t)/saveInd),1);

% Assign initial values
zR_S(:,1) = zR;
H_S(:,1)  = H;
Q         = zeros(length(x)+1,1);
k         = k*ones(1,length(x));

HW = H - 75;
HW = HW.*(HW >= 0);

%% Main loop
% Initialize counter for saving output
jj = 1;

% Loop
for ii = 1:length(t)
    % Compute glacier surface elevation 
    z = zR + H;
        
    % Compute current ELA
    ELA = 2500;
    
    % Evaluate net accumulation/ablation 
    b = gamma*(z - ELA);
    
    % Cap accumulation
    b(b > bmax) = bmax;
    
    % Compute ice surface slope
    dzdx = diff(z)/dx;
    
    % Get box-centered heights
    Hm = (H(1:end-1) + H(2:end))/2;
    HWm = (HW(1:end-1) + HW(2:end))/2;
    HWm = Hm.*(HWm >= Hm) + HWm.*(HWm < Hm);
    
    % Compute base slope
    dzRdx = diff(zR)/dx;
    
    % Compute basal shear stress [Pa]
    tauB = rhoI*g*Hm.*dzRdx;
        
    % Compute sliding speed  
    Ne = rhoI*g*Hm - rhoW*g*HWm;
    usl = zeros(size(Hm));
    usl(Hm > 0) = (usl0*tauB(Hm > 0).^2)./Ne(Hm > 0);
    
    % Compute flux
    Q(2:end-1) = usl.*Hm - A*((rhoI*g*dzdx).^3).*((Hm.^5)/5);
    
    % Error trap (stops simulation if there is a numeric crash)
    if any(isnan(Q))
        error('NaN in flux array Q')
    end
    
    % Compute flux gradient 
    dQdx = diff(Q)/dx; 
    
    % Compute thickness rate 
    dHdt = b - dQdx';
    
    % Compute bulk bed erosion rate
    eDot = edot0 * usl;

    % Update glacier thickness
    H = H + dHdt*dt;
    
    % Remove negative thickness
    H = H.*(H >= 0);
    
    % Determine melt rate
    m = 2*sin((2*pi/1)*t(ii))*0.06*(3.5 - z/1000).*(H > 0);
    m = -m.*(m < 0); % only negative change in height is melt and fix sign
    m = (rhoI/rhoW)*m; % scale volume by relative density of ice/water 
    m = m*(1/p); % height scaled by porosity
    mm = (m(1:end-1) + m(2:end))/2;
    
    % Water table elevation
    zW = zR + HW; 
    zW = z.*(zW >= z) + zW.*(zW < z); % limit to glacier surface height
    
    % Head gradient
    dzWdx = diff(zW)/dx;
    
    % Flux law for water 
    km = (k(1:end-1) + k(2:end))/2;
    uW = -km.*dzWdx .*sqrt(g*(HWm*p));
    Q_W = uW.*(HWm*p);
    
    % Fix end conditions for water flux
    Q_W = [0 Q_W].*(H > 0);
    
    % Check for freezing of channels
%     fc = sin(2*pi*t(ii));
%     if fc < 0
%         Q_W = zeros(size(Q_W));
%     end
    
    % Flux gradient
    dQWdx = diff(Q_W)/dx;
    
    % Change in water thickness
    dHWdt = dQWdx + mm;
    
    % Update water cells
    HW = [HW(1:end-1) + dHWdt*dt HW(end) + dHWdt(end)*dt];
    HW = HW.*(HW > 0); % remove negative heights
    
    % Restrict water height to glacier thickness
    HW = H.*(HW >= H) + HW.*(HW < H);
    
    % Determine change in hydraulic conductivity
    dk_flow  =  scl_flow  * Q_W;
    dk_flow  =  dk_flow.*(dk_flow >= 0);
    dk_slide =  scl_slide * [0 usl];
    dk_close = -scl_close * H.^3;
    dk = dk_flow + dk_slide + dk_close;
    dk = dk.*(H > 0);
    
    % Update hydraulic conductivity
    k = k + dk;
    k = k.*(k >= 1e-8) + 1e-8*(k < 1e-8); % enforces lower limit
    k = k.*(k <= 1e2)  + 1e2 *(k > 1e2) ; % enforces upper limit
    
    % Update valley 
    zR = zR - [0 eDot*dt];
    
    if ii > 10000
        a = 1;
    end
    
    % Check if save point and log data
    if mod(ii,saveInd) == 1
        % Compute ice volume and save
        V(jj) = trapz(x,H);
        Vw(jj) = trapz(x,HW);

        % Save output 
        ELA_S(jj)   = ELA;
        zR_S(:,jj)  = zR;   % rock elevation [m]
        H_S(:,jj)   = H;    % glacier thickness [m]
        HW_S(:,jj)  = HW;
        Q_S(:,jj)   = Q';   % discharge [m^3/yr]
        QW_S(:,jj)  = Q_W'; % water discharge [m^3/yr]
        usl_S(:,jj) = usl';
        k_S(:,jj)   = k';
        m_S(:,jj)   = m';
        hydro(jj)   = Q_W(find(H > 0,1,'last')-1);
        
        % Increment counter
        jj = jj + 1;
    end
    
    % Update progress bar
    if mod(ii,100) == 1
        progressbar(ii/length(t))
    end
end

% Clean up progress bar
progressbar(1)

%% Plots
% Animation
if 1 % set to 1 to enable animation or 0 to disable
figure
set(gcf,'Position',[100 100 1200 600])
M = [];
for ii = 1:floor(length(t)/saveInd)
    % Find ELA position
    ELAind = find(zR_S(:,ii)+H_S(:,ii) < ELA_S(ii),1,'first');
    ELApos = x(ELAind)/1000;
    
    % Plot glacier
    subplot(3,2,[1,2])
    title('Glacier','Fontsize',14)
    % Bed
    fill([zeros(size(x)) x/1000],...
        [zeros(size(x)) zR_S(:,ii)'],'k')
    hold on
    % Glacier
    fill([x/1000 fliplr(x/1000)],...
        [zR_S(:,ii)' fliplr(zR_S(:,ii)'+H_S(:,ii)')],'c')
    % Water table
    fill([x/1000 fliplr(x/1000)],...
        [zR_S(:,ii)' fliplr(zR_S(:,ii)'+HW_S(:,ii)')],'b')
    hold on
    % ELA elevation
    plot([min(x/1000) max(x/1000)],[ELA_S(ii) ELA_S(ii)],'--b')
    % ELA position marker
    plot([ELApos ELApos],[0 4200],'--b')
    hold off
    % Write current ELA on plot
    tH = text('String',['ELA: ' num2str(ELA_S(ii),4) ' m']);
    tH.Position = [150 3750];
    tH.Color = 'b';
    % Write current date on plot
    tH = text('String',['Time: ' num2str(t(ii)*saveInd,4) ' yr']);
    tH.Position = [150 500];
    % Set limits and axes
    ylim([0 4200])
    ylabel('Elevation [m]')
    xlabel('Position [km]')
    
    % Plot flux
    subplot(3,2,3)
    % Flux
    plot(x/1000,Q_S(2:end,ii))
    hold on
    % ELA position marker
    plot([ELApos ELApos],[0 1e5],'--b')
    hold off
    %ylim([0 1e5])
    ylabel('Flux [m^3/yr]')
    xlabel('Position [km]')
    
    % Sliding speed
    subplot(3,2,4)
    plot(x/1000,[0; usl_S(:,ii)])
    ylabel('Sliding [m/yr]')
    xlabel('Position [km]')
    
    % Conductivity
    subplot(3,2,5)
    plot(x/1000,k_S(:,ii))
    ylabel('Conductivity [m/yr]')
    xlabel('Position [km]')
    %ylim([0 5e-3])
    
    % Water flux
    subplot(3,2,6)
    plot(x/1000,QW_S(:,ii))
    ylabel('Flux [m^3/yr]')
    xlabel('Position [km]')
    %ylim([0 5e-5])
    
    % Save frame
%     M = [M getframe(gcf)];
    pause(0.01)    
end
end

% Make movie
% v = VideoWriter('glacier.m4v','MPEG-4');
% open(v)
% for ii = 1:length(M)
%     writeVideo(v,M(ii))
% end
% close(v)

% Peak sliding speed
figure
plot(t(1:saveInd:end),max(usl_S))
ylabel('Peak Sliding Speed [m/yr]')
xlabel('Time [yr]')
title('Sliding Speed')

% Change in water table
figure
plot(x/1000,HW_S(:,end) - HW_S(:,1))
ylabel('\Delta H [m]')
xlabel('Position [km]')
title('Change in Water Table')

% Conductivity
kmin = zeros(length(x),1);
kmax = zeros(length(x),1);
for ii = 1:length(x)
    kmin(ii) = min(k_S(ii,:));
    kmax(ii) = max(k_S(ii,:));
end
figure
hold on
plot(x/1000,kmin)
plot(x/1000,kmax)
hold off
legend('Min','Max')
ylabel('Conductivity [m/yr]')
xlabel('Position [km]')
title('Min/Max Conductivity')

% Hydrograph
figure
plot(t(1:saveInd:end),hydro)
xlabel('Time [yr]')
ylabel('Flow [m^3/yr]')
title('Hydrograph')

% Volume of water
figure
plot(t(1:saveInd:end),Vw)
xlabel('Time [yr]')
ylabel('Volume [m^3]')
title('Water Volume')
