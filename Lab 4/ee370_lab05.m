% Read the s1p files
fns = 'shorted_cable.s1p';
[f(:),S11_real(:),S11_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11_short = complex(S11_real, S11_imag);

fns = 'opened_cable.s1p';
[f(:),S12_real(:),S12_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11_open = complex(S12_real, S12_imag);
%%
c0 = 3 * 10^8; 
N = 101;
fmax = 1.5E9; 

% Testing for cable
length = 0.21; 
Z0 = 50; 
VF = 0.7;
alpha = 0; 

% Frequency definitions
nfreq = 2*N + 1; % includes negative frequencies and zero
fstrt = -fmax;
fstop = fmax;
fstep = (fstop - fstrt) / (2*N);
freq = fstrt:fstep:fstop; 
omega = 2 * pi * freq;
omega_step = 2*pi*fstep;

% Time definitions
ntime = 1515; %Addjusting plot resolution
tstrt = -1E-9;
tstop = 4E-9;
tstep = (tstop - tstrt) / (ntime - 1);
time = tstrt:tstep:tstop;

W0 = N; 

f_t_open = zeros(1, ntime);
f_t_short = zeros(1, ntime);

% Define Gaussian window parameters
epsilon = 0.11512925465;
p = 30;

for i = 1:ntime
    t = time(i);
    sum_value_open = 0;
    sum_value_short = 0;
    
    for n = -N:N
        W = exp(epsilon * -p * (n/N).^2);

        % S11 Open
        S11_open_windowed = S11_open .* W;
        sum_value_open = sum_value_open + S11_open_windowed(n + N + 1) * exp(1j * n * (2 * pi * omega_step) * t);

        % S11 Short
        S11_short_windowed = S11_short .* W;
        sum_value_short = sum_value_short + S11_short_windowed(n + N + 1) * exp(1j * n * (2 * pi * omega_step) * t);
    end
    
    f_t_open(i) = (1 / W0) * sum_value_open;
    f_t_short(i) = (1 / W0) * sum_value_short;
end

S11_open = f_t_open;
S11_short = f_t_short;

% Plotting Pulse Response
figure;
plot(time*1E9, real(S11_open));
xlabel('Time (ns)');
xlim([-1, 4])
xticks(-1:1:4)
yticks(-1:0.5:1)
ylim([-1,1])
ylabel('Pulse Amplitude (linear scale)');
title('S11 Derived Pulse Response (Open)');

% For S11 OPen
S11_mag = mag2db(abs(S11_open));
S11_phase = rad2deg(angle(S11_open));

figure;
title('S11 Derived Pulse Response (Open)');
yyaxis left
plot(time*1E9, S11_mag); 
xlim([-1, 4])
xticks(-1:1:4)

ylim([-70,10])
yticks(-70:10:10)
xlabel('Time (ns)');
ylabel('Pulse Amplitude (in dB)');

hold on
yyaxis right
plot(time*1E9, S11_phase);
ylim([-0,150])
yticks(0:50:150)
xlabel('Time (ns)');
ylabel('Pulse Phase (in Degrees)');

% For S11 Short
figure;
plot(time*1E9, real(S11_short));
xlabel('Time (ns)');
xlim([-1, 4])
xticks(-1:1:4)
ylim([-1,1])
yticks(-1:0.5:1)
ylabel('Pulse Amplitude');
title('S11 Derived Pulse Response (Short)');

S11_mag = mag2db(abs(S11_short));
S11_phase = rad2deg(angle(S11_short));

figure;
title('S11 Derived Pulse Response (Short)');
yyaxis left
plot(time*1E9, S11_mag); 
xlim([-1, 4])
xticks(-1:1:4)

ylim([-70,10])
yticks(-70:10:10)
xlabel('Time (ns)');
ylabel('Pulse Amplitude (in dB)');

hold on
yyaxis right
plot(time*1E9, S11_phase);
ylim([-0,150])
yticks(0:50:150)
xlabel('Time (ns)');
ylabel('Pulse Phase (in Degrees)');