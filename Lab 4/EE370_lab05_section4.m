%% 4.1 Modify for Touchstone
c0 = 3 * 10^8; 
N = 101;
fmax = 1.5E9;

% Read data from the file
fns = 'cascade_of_cables.s1p';
[f(:),S11_real(:),S11_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11_short = complex(S11_real, S11_imag);

% Frequency definitions
nfreq = N; % includes negative frequencies and zero
fstrt = 0; % start frequency from the data
fstop = fmax; % stop frequency from the data
fstep = (fstop - fstrt) / (N-1);
omega = 2 * pi * f;
omega_step = 2*pi*fstep;

% Time definitions
ntime = 202;
tstrt = -5E-9;
tstop = 1/fstep;
tstep = (tstop - tstrt) / (ntime - 1);
time = tstrt:tstep:tstop;

% Define Gaussian window parameters
epsilon = 0.11512925465;
p = 30;
W0 = nfreq;

W = [];
for n = 1:202
    W = [W, exp(-epsilon*p*((n-nfreq)/nfreq^2))];
end

f_t_short = zeros(1, ntime);

for i = 1:ntime
    t = time(i);
    sum_value_short = 0;
    
    for n = 1:nfreq
        % S11 Short
        S11_short_windowed = S11_short .* W;
        sum_value_short = sum_value_short + S11_short_windowed(i) * exp(1j * n * (2 * pi * omega_step) * t);
    end
    
    f_t_short(i) = (1 / W0) * sum_value_short;
end

S11_short = f_t_short;

% Interpolating S11_short for a smoother curve
time_interp = linspace(tstrt, tstop, 1515);
S11_short_interp = interp1(time, S11_short, time_interp, 'spline');

% For S11 Short
figure;
plot(time_interp*1E9, real(S11_short_interp));
xlabel('Time (ns)');
xlim([-5, 135])
ylim([-1,1])
yticks(-1:0.5:1)
ylabel('Pulse Amplitude (Linear Scale)');
title('S11 Derived Pulse Response (Short)');

%S11_mag_s_interp = mag2db(abs(S11_short_interp));
S11s_mag = mag2db(sqrt((real(S11_short_interp)).^2 + (imag(S11_short_interp )).^2)); %mag2db(abs(f_t));
S11s_phase = rad2deg(atan2(imag(real(S11_short_interp)), real(real(S11_short_interp)))); %rad2deg(unwrap(angle(f_t)));
%S11_phase_s_interp = rad2deg(unwrap(angle(S11_short_interp)));
%S11_phase_s_unwrapped = rad2deg((angle(S11_phase_s_interp)));

figure;
title('Interpolated S11 Derived Pulse Response (Short)');
yyaxis left
plot(time_interp*1E9, S11s_mag);
xlim([-1, 4]);
xticks(-1:1:4)
xlabel('Time (ns)');
ylabel('Pulse Amplitude (in dB)');

hold on
yyaxis right
plot(time_interp*1E9, S11s_phase);
ylim([0,150])
yticks(0:50:150)
xlabel('Time (ns)');
ylabel('Pulse Phase (in Degrees)');
