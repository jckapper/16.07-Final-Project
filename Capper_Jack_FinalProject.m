clear all;
close all;

% SET ALL CONSTANTS
AU = 149.598e9; %m/AU
day = 60*60*24; %s/day
G = 6.67430e-11*(day^2)/(AU^3); %AU^3/kg*day^2
M_sun = 1.98847e30; %kg
mu = G*M_sun; %AU^3/day^2;
a_earth = 1;  %AU
r_earth = 4.26e-05; %AU
v_earth = 1.720e-02; %AU/day
a_mars  = 1.523679; %AU
r_mars  = 2.27e-05;  %AU
v_mars  = 1.394e-02; %AU/day

% -----------
% PROBLEM 1 |
% -----------

% INITIALIZE ODE FUNCTION AND OPTIONS
orbital_pos = @(t,X) [X(3); X(4);-mu*X(1)/((X(1)^2+X(2)^2)^1.5);-mu*X(2)/((X(1)^2+X(2)^2)^1.5)];
options = odeset('RelTol',1e-6, 'Abstol',1e-6);        

% EARTH VARIABLES
r_0_earth = [a_earth; 0];% Earth position at perhelion (Planar so two dimensional)
v_0_earth = [0; v_earth];% At perigee, velocity is entirely tangential
initial_earth = [r_0_earth; v_0_earth];
tspan = [0:0.1:2*365]; 
[t_earth, X_earth] = ode45(orbital_pos, tspan, initial_earth, options);

% MARS VARIABLES
r_0_mars = [a_mars; 0];% Mars position at perhelion (Planar so two dimensional)
v_0_mars = [0; v_mars];% At perigee, velocity is entirely tangential
initial_mars = [r_0_mars; v_0_mars];
tspan = [0:0.1:3*365]; 
[t_mars, X_mars ] = ode45(orbital_pos, tspan, initial_mars , options);

% EARTH AND MARS ORBITS PLOT
figure
hold on
title('Hohmann Transfer from Earth to Mars')
xlabel('x-position (AU)')
ylabel('y-position (AU)')
xlim([-2.5 2.5])
ylim([-2 2])
plot(X_earth(:,1),X_earth(:,2),'linewidth',2)
plot(X_mars(:,1),X_mars(:,2),'linewidth',2)
plot(0,0,'.y', 'MarkerSize', 25)% Sun not to scale, simply for illustration
legend({'Earth','Mars','Sun'},'Location','southwest')
hold off

% PERIOD OF EARTH CALCULATION
x_e1 = X_earth(:,1);
pks = findpeaks(x_e1);
pk_idx = find(x_e1==pks(1));
period_e_obs = t_earth(pk_idx(1))% Period in days
period_e_thry = sqrt(4*pi^2*r_0_earth(1)^3/mu)% Period in days
error_e = abs((period_e_thry-period_e_obs)/period_e_thry)
T_earth = period_e_obs;

% PERIOD OF MARS CALCULATION
x_m1 = X_mars(:,1);
pks = findpeaks(x_m1);
pk_idx = find(x_m1==pks(1));
period_m_obs = t_mars(pk_idx(1))% Period in days
period_m_thry = sqrt(4*pi^2*r_0_mars(1)^3/mu)% Period in days
error_m = abs((period_m_thry-period_m_obs)/period_m_thry)
T_mars = period_m_obs;

% --------------------------------
% PROBLEM 2 INITIAL CALCULATIONS |
% --------------------------------

% HOHMANN TRANSFER VARIABLES
conversion_kms = (AU/day)/1000; %AU/Day to km/s
v_pi = sqrt(mu*(2/a_earth - 2/(a_earth+a_mars)));
v_alpha  = sqrt(mu*(2/a_mars -  2/(a_earth+a_mars)));
v_alpha_kms = v_alpha*conversion_kms;
v_pi_kms = v_pi*conversion_kms;

% DELTA_V REQUIRED TO TRANSFER
delta_V_perogee = v_pi - v_earth;
delta_V_apogee  = v_mars - v_alpha;

% HOHMANN TRANSFER PERIOD
a_hohmann = (a_earth + a_mars)/2;
T_hohmann = 2*pi*sqrt((a_hohmann^3)/mu);
T_transfer = T_hohmann/2; % Period in days

% CALCULATE HOHMANN TRANSFER
r_0_hohmann = [a_earth; 0];
v_0_hohmann = [0; v_pi];
initial_h = [r_0_hohmann; v_0_hohmann];
[t_hohmann,X_hohmann] = ode45(orbital_pos, tspan, initial_h, options);

% PLOT TRANSFER ORBIT
figure
hold on
title('Hohmann Transfer from Earth to Mars')
xlabel('x-position (AU)')
ylabel('y-position (AU)')
xlim([-2.5 2.5])
ylim([-2 2])
plot(X_earth(:,1),X_earth(:,2),'linewidth',2)
plot(X_mars(:,1),X_mars(:,2),'linewidth',2)
plot(X_hohmann(:,1), X_hohmann(:,2),'--','linewidth',1.5)
plot(0,0,'.y', 'MarkerSize', 25)% Sun not to scale, simply for illustration
legend({'Earth','Mars','Hohmann Transfer','Sun'},'Location','southwest')
hold off

%----------------
% PROBLEM 2a/2b |
%----------------

theta_earth_i = 104;  %degrees
theta_mars_i  = 216;  %degrees

% ANGLE REQUIRED FOR SUCCESSFUL TRANSFER ORBIT
omega_earth = 360/T_earth;
omega_mars  = 360/T_mars;
theta_req_em = 180-omega_mars*T_transfer;
theta_req_me = 180-omega_earth*T_transfer;

% LOCATION OF PLANET ALONG ORBIT AT TIME t
t = 0:1:6*365;
theta_earth = t*omega_earth + theta_earth_i;
theta_mars  = t*omega_mars  + theta_mars_i ;

% SUCCESSFUL LAUNCH TIMES FROM EARTH CALCULATION AND PLOT
figure
hold on
dept_times_earth = [];
arrive_times_mars = [];
for n=-4:1:0
    t = ((-67.5497+(360*n))/(omega_mars-omega_earth));
    dept = 2020+(t/365.25);% Time converted to years
    dept_times_earth = [dept_times_earth dept];
    arrival = dept+T_transfer/365;
    arrive_times_mars = [arrive_times_mars arrival];
    scatter(dept,arrival, 'filled')
end
xticks(2020:1:2030)
xlabel('Earth Departure Time (Years)')
ylabel('Mars Arrival Time (Years)')
title('Earth to Mars Transit Dates')
hold off
launch_time_interval = dept_times_earth(1)-dept_times_earth(2);
%-------------
% PROBLEM 2c |
%-------------

% SUCCESSFUL LAUNCH TIMES FROM MARS CALCULATION AND PLOT
figure
hold on
dept_times_mars = [];
arrive_times_earth = [];
for n=1:1:4
    t = ((38.7538+(360*n))/(omega_earth-omega_mars));
    dept = 2020+(t/365.25);% Time converted to years
    dept_times_mars = [dept_times_mars dept];
    arrival = dept+T_transfer/365;
    arrive_times_earth = [arrive_times_earth arrival];
    scatter(dept,arrival,'filled')  
end
xticks(2020:1:2030)
xlabel('Mars Departure Time (Years)')
ylabel('Earth Arrival Time (Years)')
title('Mars to Earth Transit Dates')
hold off

%-------------------
% PROBLEM 2d/2e/2f |
%-------------------

% REQUIRED STAY ON MARS
t_on_mars = dept_times_mars(1)-arrive_times_mars(length(arrive_times_mars))

% MINIMUM ROUND TRIP TIME
t_round_trip = arrive_times_earth(1)-dept_times_earth(length(dept_times_earth))

% DELTA V's OVER ENTIRE TRIP
conversion = AU/day;    %AU/Day to m/s
delta_V_launch = delta_V_perogee*conversion;
delta_V_arrival = delta_V_apogee*conversion;
delta_V_return_launch = -delta_V_apogee*conversion;
delta_V_return_arrival = -delta_V_perogee*conversion;
delta_V_total = 2*abs(delta_V_launch) + 2*abs(delta_V_arrival);

%----------------
% PROBLEM 3a/3b |
%----------------

% LAMBERT SPLASH PLOT INITIALIZATION
figure
xline(0,'-','HandleVisibility','off');
yline(0,'-','HandleVisibility','off');
hold on
plot(X_earth(:,1),X_earth(:,2),'DisplayName','Earth','LineWidth',2)
title('Earth to Mars Lambert Splash Transfers')
xlabel('X Position (AU)')
ylabel('Y Position (AU)')
plot(X_mars(:,1),X_mars(:,2),'DisplayName','Mars','LineWidth',2)
plot(0,0,'.y', 'DisplayName', 'Sun', 'MarkerSize', 25)% Sun not to scale, simply for illustration
lgd = legend('Location','northeast');
lgd.FontSize = 4;

% INITIALIZE VALUES
deltas = [-68.8,-46,-23,0,23,46,68.8];
final_values = [];
delta_V_lambert = 2*delta_V_perogee;
for delta=deltas
    
    % SET TSPAN TO GET CORREACT PLOT
    if delta == 0
        tspan = [0:0.1:2.6*365];
    elseif delta == -23 || delta == 23 
        tspan = [0:0.1:2.2*365];
    else
        tspan = [0:0.1:2*365];
    end
    
    % PLOT LAMBERT SPLASH ORBIT
    r_0_lambert = [a_earth; 0];% Mars position at perigee (Planar so two dimensional)
    v_0_lambert = [delta_V_lambert*sind(delta); v_earth+(delta_V_lambert*cosd(delta))];% DeltaV is twice Hohmann and initial condition varys based on the angle delta
    initial_h = [r_0_lambert; v_0_lambert];
    [t_lambert,X_lambert] = ode45(@(t,x) orbital_pos(mu,x), tspan, initial_h, options);
    plot(X_lambert(:,1),X_lambert(:,2),'DisplayName',sprintf('Delta=%f', delta))
    
    % EXTRACT POSITION VALUES
    x_vals = X_lambert(:,1);
    y_vals = X_lambert(:,2);
    
    % FIND INTERSECTION POINT 
    min_idx = find(abs(sqrt(X_lambert(:,1).^2 + X_lambert(:,2).^2)-a_mars) <= 1e-3, 1);
    max_idx = find(abs(sqrt(X_lambert(:,1).^2 + X_lambert(:,2).^2)-a_mars) <= 1e-3, 1, 'last');
    
    % CALCULATE THETA TRAVELED, TOF, AND THETA REQUIRED
    if delta == 68.8 || delta == -68.8
        if delta == -68.8
            intersect_1 =  X_lambert(max_idx,1:2);
        else
            intersect_1 =  X_lambert(min_idx,1:2);
        end
        plot([0 intersect_1(1)],[0 intersect_1(2)],'--b','HandleVisibility','off')
        TOF_1 = t_lambert(min_idx);
        theta_trvld_1 = atand(intersect_1(2)/intersect_1(1))+180;
        theta_offset_1 = mod(theta_trvld_1 - (360/T_mars)*TOF_1, 360);
        final_values = [final_values; [delta, TOF_1, theta_trvld_1, theta_offset_1]];
    else
        intersect_1 = X_lambert(min_idx,1:2);
        if intersect_1(1)>0
            theta_trvld_1 = atand(intersect_1(2)/intersect_1(1));
            if theta_trvld_1<0
                theta_trvld_1 = theta_trvld_1+360;         
            end
        else
            theta_trvld_1 = atand(intersect_1(2)/intersect_1(1))+180;
        end
        
        intersect_2 = X_lambert(max_idx,1:2);
        if intersect_2(1)>0
            theta_trvld_2 = atand(intersect_2(2)/intersect_2(1));
            if theta_trvld_2<0
                theta_trvld_2 = theta_trvld_2+360;               
            end
        else
            theta_trvld_2 = atand(intersect_2(2)/intersect_2(1))+180;
        end
        
        TOF_1 = t_lambert(min_idx);
        TOF_2 = t_lambert(max_idx);
        theta_offset_1 = mod(theta_trvld_1 - (360/T_mars)*TOF_1, 360);
        theta_offset_2 = mod(theta_trvld_2 - (360/T_mars)*TOF_2, 360);
        final_values = [final_values; [delta, TOF_1, theta_trvld_1,theta_offset_1]];
        final_values = [final_values; [delta, TOF_2, theta_trvld_2,theta_offset_2]];
        plot([0 intersect_1(1)],[0 intersect_1(2)],'--b','HandleVisibility','off')
        plot([0 intersect_2(1)],[0 intersect_2(2)],'--b','HandleVisibility','off')
    end
end    
hold off

%----------------
% PROBLEM 3c/3d |
%----------------

% CALCULATE THETA REQUIRED FOR LAMBERT SPLASH ORBITS
launch_info = [];
t_span = 700;
for i = 1:size(final_values)
    angle_offset = final_values(i,4);
       TOF = final_values(i,2);
    for j = 1:1:t_span
        current_date = 2020 + j/T_earth;
        current_angle = theta_mars(j)-theta_earth(j);
        if abs(mod((current_angle-angle_offset),360)) < 0.5
            launch_info = [launch_info; [current_date, TOF]];
        end
    end
end

% PORKCHOP PLOT
figure
title('Launch Dates and Flight Durations');
xlabel('Launch Date (Year)');
ylabel('Transit Duration (Days)');
hold on
scatter(launch_info(:,1),launch_info(:,2),'filled');
scatter(dept_times_earth(length(dept_times_earth)), T_transfer, 100, 'p')
hold off