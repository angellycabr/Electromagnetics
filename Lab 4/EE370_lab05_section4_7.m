% Read the s1p files
fns = 'cascade_of_cables.s1p';
[f(:),S11_real(:),S11_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11_short = complex(S11_real, S11_imag);
%% Plotting the Results

% Specify frequencies
nfreq = 202;
fstrt = 0.0*10^9;
fstop = 1.5*10^9;
fstep = (fstop-fstrt)/(nfreq-1);
freq = -fstop:fstep:fstop;
omega = 2*pi*freq;
omega_step = 2*pi*fstep;

% Specify times
ntime = 1001;
tstrt = -5 * 10^-9;
tstop = 1/fstep; %65*10^-9;
tstep = (tstop-tstrt)/(ntime-1);
time = tstrt:tstep:tstop;

u_t = zeros(size(time));
u_t(time*1E9 >= 2) = 1;

% Compute gamma values
c0 = 3*1E8;
VF = 0.7;  
alpha = 0; 

epsilon = 0.11512925465;
p = 30;

% Compute f(t) using SIFT for S11_open(n)
f_t = zeros(1, ntime);
for k = 1:ntime
    sum_val = 0;

    for n = 1:nfreq
        sum_val = sum_val + S11_short(n) * exp(epsilon * -p * (n/nfreq).^2) * exp(1j * n * omega_step * time(k));
    end
    conv_result = conv(f_t, u_t, 'same');
    f_t(k) = (1/nfreq) * sum_val;
end

% Plot the pulse response
figure;
plot(time*1E9, real(f_t));
xlim([tstrt*1E9, (tstop*1E9)+1])
xticks(tstrt*1E9:10:(tstop*1E9)+1)
%yticks(-0.5:0.5:0.5)
ylim([-0.5,0.5])
xlabel('Time (s)');
ylabel('Pulse Response (Linear Scale)');
title('S11 Derived Pulse Response');

% Finding magntiude and phase
S11p_mag = mag2db(sqrt((real(f_t)).^2 + (imag(f_t)).^2)); %mag2db(abs(f_t));
S11p_phase = rad2deg(atan2(imag(real(f_t)), real(real(f_t)))); %rad2deg(unwrap(angle(f_t)));

figure;
title('S11 Derived Pulse Response');
yyaxis left
plot(time*1E9, S11p_mag); 
xlim([-5, 10])
xticks(-5:5:10)
ylim([-60,10])
yticks(-60:10:10)
xlabel('Time (ns)');
ylabel('Pulse Amplitude (in dB)');

hold on
yyaxis right
plot(time*1E9, S11p_phase);
xlabel('Time (ns)');
ylim([0,200])
yticks(0:50:200)
ylabel('Pulse Phase (in Degrees)');
%%
%u_t = zeros(size(time));
%u_t(time*1E9 >= 2) = 1;

% Perform convolution between f_t and Heaviside function
%conv_result = conv(real(f_t), u_t, 'same');

% Plot the convolved result over [0, 5] ns time interval
figure;
plot(time*1E9, real(conv_result));
xlim([0, 5]) % Limit the x-axis to the [0, 5] ns range
xticks(0:1:5)
ylim([-1,-0.5])
yticks(-1:0.1:-0.5)
xlabel('Time (ns)');
ylabel('Convolution Result');
title('Convolution of f_t and Heaviside Function');
%%
% Define the rectangular pulse function over the time range
half_width = 0.4 * 10^-9 / 2; % Half width of the pulse
rect_pulse = zeros(size(time));
rect_pulse((time >= -half_width) & (time <= half_width)) = 1;

% Perform convolution between f_t and rectangular pulse
conv_result = conv(f_t, rect_pulse, 'same');

% Plot the convolved result over [0, 5] ns time interval
figure;
plot(time*1E9, real(conv_result));
xlim([0, 5]) % Limit the x-axis to the [0, 5] ns range
xticks(0:1:5)
xlabel('Time (ns)');
ylabel('Convolution Result');
title('Convolution of f_t and Rectangular Pulse');