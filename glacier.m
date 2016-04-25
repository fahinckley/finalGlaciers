%--------------------------------
% Simulation of an eroding valley under a glacier
%--------------------------------
% Franklin Hinckley
% 6 April 2016
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
k = 1000/315569.26; % initial hydraulic conductivity [m/yr]
p = 0.05; % porosity of glacier
scl_flow = 1e-5; % scale factor for conductivity as a function of water flow
scl_slide = 1e-16; % scale factor for conductivity as a function of sliding
scl_close = 1e-18; % scale factor for conductivity as a function of ice thickness

%% Set up initial valley shape
dx = 1000; % [m]
len = 200*1000; % [m]
x = 0:dx:len;

SR = 20/1000;  % [m/km]

zMax = 3000; % [m]
zR = zMax - SR*x;

%% Time array
dt   = 1/365/24; % [yr]
tSim = 5;    % [yr]
t = 0:dt:tSim;

P_ELA = 1; % period for variations in ELA [yr]
P_W   = 1;   % period for variations in water table [yr]

%% Initial glacier shape
load Hsteady
H = Hsteady'; % steady state glacier from original code

%% Allocate output
% Specify how often to save
saveInd = 100;

% Allocate output
zR_S  = zeros(length(x)  , floor(length(t)/saveInd));
H_S   = zeros(length(x)  , floor(length(t)/saveInd));
HW_S  = zeros(length(x)  , floor(length(t)/saveInd));
k_S   = zeros(length(x)  , floor(length(t)/saveInd));
usl_S = zeros(length(x)-1, floor(length(t)/saveInd));
Q_S   = zeros(length(x)+1, floor(length(t)/saveInd));
QW_S  = zeros(length(x)+1, floor(length(t)/saveInd));
m_S  = zeros(length(x), floor(length(t)/saveInd));
ELA_S = zeros(floor(length(t)/saveInd),1);
V     = zeros(floor(length(t)/saveInd),1);

% Assign initial values
zR_S(:,1) = zR;
H_S(:,1)  = H;
Q         = zeros(length(x)+1,1);
k         = k*ones(1,length(x));
% HW        = H - 175;
%HW = zeros(size(H));
HW = 100*ones(size(H));

HW = HW.*(HW >= 0);

%% Main loop
% Initialize counter for saving output
jj = 1;

% Loop
for ii = 1:length(t)
    % Compute glacier surface elevation 
    z = zR + H;
        
    % Compute current ELA
    ELA = 2500 - 1500*sin((2*pi/P_ELA)*t(ii));
    
    % Evaluate net accumulation/ablation 
    b = gamma*(z - ELA);
    
    % Cap accumulation
    b(b > bmax) = bmax;
    
    % Compute ice surface slope
    dzdx = diff(z)/dx;
    
    % Get box-centered heights
    Hm = (H(1:end-1) + H(2:end))/2;
    HWm = (HW(1:end-1) + HW(2:end))/2;
    
    % Compute base slope
    dzRdx = diff(zR)/dx;
    
    % Compute basal shear stress [Pa]
    tauB = rhoI*g*Hm.*dzRdx;
    
    % Compute water table level
    %Hwtl = 65 + 10*sin((2*pi/P_W)*t(ii));
    %Hwtl = 75;
    %Hw = Hm - Hwtl;
    %Hw = Hw.*(Hw > 0);
    
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
    m = b.*(H > 0)*dt; % assumes no sublimation
    m = -m.*(m < 0); % only negative change in height is melt and fix sign
    m = (rhoI/rhoW)*m; % scale volume by relative density of ice/water 
    m = m*(1/p); % height scaled by porosity
    
    % Water table elevation
    zW = zR + HW; 
    zW = z.*(zW >= z) + zW.*(zW < z);
    
    % Head gradient
    dzWdx = diff(zW)/dx;
    
    % Flux law for water 
    km = (k(1:end-1) + k(2:end))/2;
    Q_W = -km.*dzWdx;
    
    % Fix end conditions for water flux
    Q_W = [0 Q_W Q_W(end)];
    
    % Freeze sections above ELA
    %Q_W = Q_W.*(z <= ELA);
    
    % Flux gradient
    dQWdx = diff(Q_W)/dx;
    
    % Change in water thickness
    dHWdt = dQWdx + m;
    
    % Update water cells
    HW = HW + dHWdt*dt;
    HW = HW.*(HW > 0); % remove negative heights
    
    % Restrict water height to glacier thickness
    HW = H.*(HW >= H) + HW.*(HW < H);
    
    % Determine change in hydraulic conductivity
    dk_flow  =  scl_flow  * Q_W(1:end-1);
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
    ylim([0 1e5])
    ylabel('Flux [m^3/yr]')
    xlabel('Position [km]')
    
    % Sliding speed
    subplot(3,2,4)
    plot(x/1000,[0; usl_S(:,ii)])
    ylabel('Sliding Speed [m/yr]')
    xlabel('Position [km]')
    
    % Conductivity
    subplot(3,2,5)
    plot(x/1000,k_S(:,ii))
    ylabel('Conductivity [m/s]')
    xlabel('Position [km]')
    %ylim([0 5e-3])
    
    % Water flux
    subplot(3,2,6)
    plot(x/1000,QW_S(2:end,ii))
    ylabel('Flux [m^3/yr]')
    xlabel('Position [km]')
    %ylim([0 5e-5])
    
    % Save frame
    %M = [M getframe(gcf)];
    pause(0.01)    
end
end

% Make movie
% v = VideoWriter('glacierWater.m4v','MPEG-4');
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

% Volume vs time
Veq = (1-(1/exp(1)))*V(end); % 1 - (1/e) of equilibrium volume
eqInd = find(V > Veq,1,'first');
teq = t(eqInd); % time to Veq
figure
hold on
plot(t(1:saveInd:end),V)
%plot([t(1) t(end)],[Veq Veq],'--k')
%plot([teq teq],[0 10e7],'--k')
hold off
xlabel('Time [yr]')
ylabel('Ice Volume [m^3]')
ylim([0 1e8])
%tH = text('String',['t_{eq}: ' num2str(teq,4) ' yr']);
%tH.Position = [teq+50 1e7];
