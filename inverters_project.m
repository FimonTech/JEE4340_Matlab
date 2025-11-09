%% Inverter Project
% Waveform Parameters
Vs = 220;
F = 60;
w = 2*pi*F;
T = 1/F;
Vm = Vs*sqrt(2);

% Inverter Parameters
R = 6.5;
L = 51e-3; % induction
C = 138e-6;
Xl0 = 1i*w*L;
Xc0 = 1i/(w*C);

% FFT parameters
Fs = 2^17;          % using power of 2 cause its fast
Ts = 1/Fs;
l = round(T*Fs);
t = (0:l-1)*Ts;     % length of one period 
x1 = Vs*(2*ceil(sin(w*t))-1);
% Plotting
[output,Io] = analyze_signal(x1,l,Fs,t,16,"Square Wave inverter",1,R,L,C);
I_o = sum(Io);
figure()
plot(t,I_o,t,Io(1,:),'--'); % Plotting current output
legend(["Total Output Current","Fundamental Current"])
grid on
title("Reconstructed Load Current, I_o, R="+num2str(R)+ ...
    ",L="+num2str(L)+",C="+num2str(C))
xlabel("Time, ms")
ylabel("I_o(t), A")
        
%{
%% Appendix 1: Main code
% rectifiers_fft.m
type('inverters_project.m')

%% Appendix 2: calculation and plotting function
% construct_signal.m
type('analyze_signal.m')
%}