% Load in Data and visualize
% File Format = [time pressure flow]

global data time pressure flow 
data = importdata("tPQ_29.txt");
data = data.data;
time = data(:, 1); % units = seconds 
pressure = data(:, 2); % units = mmHg
flow = data(:, 3); %units = mL/min 

% Plot Pressure 

subplot(2, 2, 1);
plot(time, pressure, '.', 'MarkerSize', 5)
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
title('Data: Pressure vs Time')
set(gca, 'FontSize', 14)

% Plot dP/dt 

dPdt = diff(pressure) ./ diff(time);
sgfdPdt = sgolayfilt(dPdt,3,21); %smoothing the estiamtion

subplot(2, 2, 2);
plot(time(2:end), dPdt, '.-r', time(2:end), sgfdPdt, '.-b', 'MarkerSize', 10)
legend('dPdt - Raw', 'dPdt - Smoothed')
xlabel('Time (s)')
ylabel('Change in Pressure (mmHg/s)')
title('Data: Change in Pressure vs Time')
set(gca, 'FontSize', 14)

% Plot Flow

subplot(2, 2, 3);
plot(time, flow, '.', 'MarkerSize', 5)
xlabel('Time (s)')
ylabel('Flow (mL/min)')
title('Data: Flow vs Time')
set(gca, 'FontSize', 14)

% Plot dQ/dt 

dQdt = diff(flow) ./ diff(time);
sgfdPdt = sgolayfilt(dQdt,3,21); %smoothing the estiamtion

subplot(2, 2, 4);
plot(time(2:end), dQdt, '.-r', time(2:end), sgfdPdt, '.-b', 'MarkerSize', 10)
legend('dQdt - Raw', 'dQdt - Smoothed')
xlabel('Time (s)')
ylabel('Change in Flow (mL/s^2)')
title('Data: Change in Flow vs Time')
set(gca, 'FontSize', 14)

print -djpeg RawDataPlots

% Two Element Windkessel Model

Rd = 0.86;        % mmHg*(s/ml)
C  = 0.17;        % ml/mmHg
ts = [time(1) time(end)];
p0 = pressure(1);
theta0 = [Rd C];

% Pre-interpolate Q once as a function handle

Q_interp = @(t) interp1(time, flow, t, 'linear', 'extrap');

function optimize = objective_Wind2Pressure(theta, time, pressure, Q_interp, p0, ts)
    Rd = theta(1);
    C  = theta(2);
    try
        [t2, p2] = ode45(@(t,p) Wind2Pressure(t, p, Q_interp, C, Rd), ts, p0);
        Pmodint  = interp1(t2, p2, time, 'linear', 'extrap');
        optimize = sum((pressure - Pmodint).^2);
    catch
        optimize = 1e10;
    end
end

function dPdt = Wind2Pressure(t, P, Q_interp, C, Rd)
    Q    = Q_interp(t);
    dPdt = Q/C - P/(C*Rd);
end

func1 = @(theta) objective_Wind2Pressure(theta, time, pressure, Q_interp, p0, ts);
theta_optimized = fminsearch(func1, theta0);

Rd = theta_optimized(1);
C  = theta_optimized(2);

[t2, p2] = ode45(@(t,p) Wind2Pressure(t, p, Q_interp, C, Rd), ts, p0);
p2 = interp1(t2, p2, time);

% Predicting Flow 

dPdt_raw = diff(pressure) ./ diff(time);
dPdt_smooth = sgolayfilt(dPdt_raw, 3, 21);
dPdt_full = interp1(time(1:end-1) + diff(time)/2, dPdt_smooth, time, 'linear', 'extrap');

function optimize = objective_Wind2Flow(theta, pressure, dPdt_full, flow)
    Rd = theta(1);
    C  = theta(2);
    Fmodint  = pressure / Rd + C * dPdt_full;
    optimize = sum((flow - Fmodint).^2);
end

function dVdt = Wind2Flow(pressure, dPdt_full, C, Rd)
    dVdt = pressure / Rd + C * dPdt_full;
end

func1 = @(theta) objective_Wind2Flow(theta, pressure, dPdt_full, flow);
theta_optimized = fminsearch(func1, theta0);

Rd = theta_optimized(1);
C  = theta_optimized(2);
q2 = Wind2Flow(pressure, dPdt_full, C, Rd);


% Three Element Windkessel model

R = 1.2; % units = mmHg*(s/ml)
theta02 = [R C Rd];

% Predicting Pressure 

p0 = pressure(1);
ts = [time(1) time(end)];

dQdt_raw = diff(flow) ./ diff(time);
dQdt_smooth = sgolayfilt(dQdt_raw, 3, 21);
dQdt_interp = @(t) interp1(time(1:end-1), dQdt_smooth, t, 'linear', 'extrap');
Q_interp    = @(t) interp1(time, flow, t, 'linear', 'extrap');

function dPdt = Wind3Pressure(t, P, Q_interp, dQdt_interp, C, Rd, R)
    Q    = Q_interp(t);
    dQdt = dQdt_interp(t);
    dPdt = Q/C + (Q*R)/(Rd*C) - P/(C*Rd) + R*dQdt;
end

function optimize = objective_Wind3Pressure(theta, time, Q_interp, dQdt_interp, pressure, p0, ts)
    R  = theta(1);
    C  = theta(2);
    Rd = theta(3);
    try
        [t3, p3] = ode45(@(t,p) Wind3Pressure(t, p, Q_interp, dQdt_interp, C, Rd, R), ts, p0);
        Pmodint  = interp1(t3, p3, time, 'linear', 'extrap');
        optimize = sum((pressure - Pmodint).^2);
    catch
        optimize = 1e10;
    end
end

options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);
func2 = @(theta02) objective_Wind3Pressure(theta02, time, Q_interp, dQdt_interp, pressure, p0, ts);
theta_optimized3 = fminsearch(func2, theta02, options);

R  = theta_optimized3(1);
C  = theta_optimized3(2);
Rd = theta_optimized3(3);
[t3, p3] = ode45(@(t,p) Wind3Pressure(t, p, Q_interp, dQdt_interp, C, Rd, R), ts, p0);
p3 = interp1(t3, p3, time);

subplot(2, 2, 1);
plot(time, p3, '.r', time, p2, '.g', time, pressure, '.b')
legend('Pressure - Three El Model','Pressure - Two El Model','Pressure - Data')
xlabel('Time (s)'), ylabel('Pressure (mmHg)')
title('Windkessel Models versus Data')
set(gca, 'FontSize', 14)

% Pressure Error Plots 

subplot(2, 2, 2);
semilogy(time, errorWind2Pressure, time, errorWind3Pressure)
legend('Two Element', 'Three Element')
xlabel('Time (s)')
ylabel('Error in Pressure (mmHg)')
title('Windkessel Error Plot')
set(gca, 'FontSize', 14)

% Predicting Flow 

f0 = flow(1);
ts = [time(1) time(end)];

dPdt_raw = diff(pressure) ./ diff(time);
dPdt_smooth = sgolayfilt(dPdt_raw, 3, 21);
dPdt_interp = @(t) interp1(time(1:end-1), dPdt_smooth, t, 'linear', 'extrap');
P_interp    = @(t) interp1(time, pressure, t, 'linear', 'extrap');

function dQdt = Wind3Flow(t, Q, P_interp, dPdt_interp, C, Rd, R)
    P    = P_interp(t);
    dPdt = dPdt_interp(t);
    dQdt = P/(C*R*Rd) -Q/(Rd*C) - Q/(C*R) + dPdt/R;
end

function optimize = objective_Wind3Flow(theta, time, P_interp, dPdt_interp, flow, f0, ts)
    R  = theta(1);
    C  = theta(2);
    Rd = theta(3);
    try
        [t3, q3] = ode45(@(t,q) Wind3Flow(t, q, P_interp, dPdt_interp, C, Rd, R), ts, f0);
        Qmodint  = interp1(t3, q3, time, 'linear', 'extrap');
        optimize = sum((flow - Qmodint).^2);
    catch
        optimize = 1e10;
    end
end

func2 = @(theta) objective_Wind3Flow(theta, time, P_interp, dPdt_interp, flow, f0, ts);
theta_optimized3 = fminsearch(func2, theta02, options);
R  = theta_optimized3(1);
C  = theta_optimized3(2);
Rd = theta_optimized3(3);
[t3, q3] = ode45(@(t,q) Wind3Flow(t, q, P_interp, dPdt_interp, C, Rd, R), ts, f0);
q3 = interp1(t3, q3, time);

subplot(2, 2, 3);
plot(time, q3, '.r', time, q2, '.g', time, flow, '.b')
legend('Flow - Three El Model','Flow - Two El Model', 'Flow - Data')
xlabel('Time (s)'), ylabel('Flow (mL/min)')
title('Windkessel Models versus Data')
set(gca, 'FontSize', 14)

% Flow Error Plots 

subplot(2, 2, 4);
semilogy(time, errorWind2Flow, time, errorWind3Flow)
legend('Two Element', 'Three Element')
xlabel('Time (s)')
ylabel('Error in Flow (mL/min)')
title('Windkessel Error Plot')
set(gca, 'FontSize', 14)

print -djpeg ModelPlots