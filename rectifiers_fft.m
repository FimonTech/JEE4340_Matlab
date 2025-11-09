% Waveform Parameters
Vs = 120;
F = 60;
w = 2*pi*F;
T = 1/F;
Vm = Vs*sqrt(2);

% FFT parameters
Fs = 2^14;          % using power of 2 cause its fast
Ts = 1/Fs;
L = round(T*Fs);
t = (0:L-1)*Ts;     % length of one period 

%% half-wave
x1 = max(Vm*sin(w*t),0);    %function
figure(1);
construct_signal(x1,L,Fs,t,10,"Half-wave rectifier");
%%
fprintf(repmat('\n',1,8));
%% full-wave
x2 = abs(Vm*sin(w*t));    %function
figure(2);
construct_signal(x2,L,Fs,t,10,"Full-wave rectifier");

%%
fprintf(repmat('\n',1,8));
%% 3ph half-wave 
x3 = max([Vm*sin(w*t + 0*120*pi/180); ...
          Vm*sin(w*t + 1*120*pi/180);...
          Vm*sin(w*t - 1*120*pi/180)]);    %function
figure(3);
construct_signal(x3,L,Fs,t,10, "3-phase half-wave rectifier");
%%
fprintf(repmat('\n',1,12));
%% 6ph half-wave
x4 = max([  Vm*sin(w*t); ...
            Vm*sin(w*t + 1*60*pi/180); ...
            Vm*sin(w*t - 1*60*pi/180); ...
            Vm*sin(w*t + 2*60*pi/180); ...
            Vm*sin(w*t - 2*60*pi/180); ...
            Vm*sin(w*t + 3*60*pi/180)]); %function
figure(4);
construct_signal(x4,L,Fs,t,10, "6-phase half-wave rectifier");
%%
fprintf(repmat('\n',1,12));
%% 3ph full-wave 
x5 = max([Vm*sin(w*t + 0*120*pi/180); ...
          Vm*sin(w*t + 1*120*pi/180);...
          Vm*sin(w*t - 1*120*pi/180)])...
     -min([Vm*sin(w*t + 0*120*pi/180); ...
          Vm*sin(w*t + 1*120*pi/180);...
          Vm*sin(w*t - 1*120*pi/180)]); %function
figure(5);
construct_signal(x5,L,Fs,t,10,"3-phase bridge rectifier");
%%
fprintf(repmat('\n',1,12));
%% Appendix 1: Main code
% rectifiers_fft.m
type('rectifiers_fft.m')

%% Appendix 2: calculation and plotting function
% construct_signal.m
type('construct_signal.m')
