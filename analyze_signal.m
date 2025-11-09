function [output,I_o] = analyze_signal(x,l,Fs,t,N,name,p_bool,R,L,C) % updated version of construct_signal
    X = fft(x);         % this is where DSP magic happens
    P2 = abs(X/l);  
    P1 = P2(1:floor(l/2+1));
    P1(2:end-1) = 2*P1(2:end-1); % removing symmetrical frequencies at the end
    a0 = real(P1(1));  % a0 = 2 * mean(x)
    an = zeros(1,floor(l/2+1));  % Depends on how many harmonics we'll have here
    bn = zeros(1,floor(l/2+1)); 
    f = Fs*(0:floor(l/2))/l;
    [psor,~] = findpeaks(P1,f,'SortStr','none'); % Peaks calculation
    %{
    plot(f,P1) % troubleshooting plotting
    title("Single-Sided Amplitude Spectrum of X(t)")
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    figure()
    %}
    v2 = a0 * ones(size(t));
    v2_wo1 = a0 * ones(size(t)); % Waveform without the first harmonic
    %{
    disp("Harmonic number: 0"); % this is to check which harmonic worked
    disp("Harmonic ampltiude: "+num2str(a0)); % this is to check which harmonic worked
    disp("Harmonic frequency: 0")
    %}
    i = 2; % iterator for harmonics
    n = 1; % counter
    harmonics = ones(1,N);
    Harm = zeros(1,N);
    V_m = zeros(1,N);
    X_m = zeros(1,N);
    freq = zeros(1,N);
    Z = zeros(1,N);
    I_m = zeros(1,N);
    theta_n = zeros(1,N);
    I_o = zeros(N,length(t));
    while n <= N % Count to N harmonics
        an(i) = 2*real(X(i)/l);
        bn(i) = -2*imag(X(i)/l);
        if (P1(i) == psor(n)) % if point is peak
            v2 = v2 + an(i).*cos(2*pi*f(i).*t)+bn(i).*sin(2*pi*f(i).*t); % we only care about non-zero harmonics
            harmonics(n) = sqrt(an(i)^2+bn(i)^2); % Get the sum of harmonics (needed for RMS calculation)
            if (n > 1)
               v2_wo1 = v2_wo1 + an(i).*cos(2*pi*f(i).*t)+bn(i).*sin(2*pi*f(i).*t); % waveform w/o 1st harmonic
            end
            %{
            disp(an(i)); % debugging
            disp(bn(i));
            disp("Harmonic number: "+num2str(i-1)); % this shows the order of harmonic
            disp("Harmonic amplitude: "+num2str(P1(i))); % this shows the amplitude
            disp("Harmonic frequency: "+num2str(f(i))); % frequency
            %}
            if (L ~= 0 && C ~= 0 && R ~= 0)
                Harm(n) = i-1;
                V_m(n) = P1(i);
                freq(n) = f(i);
                X_m(n) = 2*pi*f(i)*L-1/(2*pi*f(i)*C); % Calculating reactance magnitude
                Z(n) = sqrt(R^2+X_m(n)^2); % Calculating impedance
                I_m(n) = P1(i)/Z(n);% Calculating current magnitude
                theta_n(n) = angle(1i*X_m(n)+R);
                I_o(n,:) = I_m(n).*sin(2*pi*f(i).*t-theta_n(n));
            n = n + 1;
            end
        end
        i = i + 1;
    end
    keys =  {'N','F(Hz)','V_m','X_L-X_C','|Z|','I_m','theta(rad)','theta(deg)'};
    output = table(Harm',freq',V_m',X_m',Z',I_m',theta_n',rad2deg(theta_n'),'VariableNames',keys);
    disp(output);
    if p_bool
        plot(t*1000,v2);
        hold on 
        plot(t*1000,v2_wo1,'--',"Color","#EDB120")
        hold off
        title("Reconstructed signal of "+name)
        xlabel("Time, ms")
        ylabel("V(t), V")
        legend(["Reconstructed V_0(t)","Reconstructed V_0(t) w/o 1st harmonic"])
        grid on
    end
    Vdc = a0;
    Vrms = sqrt(a0^2+0.5*sum(an.^2+bn.^2));
    Vac = sqrt(Vrms^2-a0^2);
    Va1rms = abs((harmonics(1)/sqrt(2)));
    
    disp("DC part:"+num2str(Vdc)+"V");
    disp("RMS:"+num2str(Vrms)+"V");
    disp("AC part:"+num2str(Vac)+"V");
    disp("Rectification Ratio:"+num2str(Vdc^2/Vrms^2));
    disp("Form Factor:"+num2str(Vrms/Vdc));
    disp("Ripple Factor:"+num2str(Vac/Vdc));
    disp("RMS Magnitude of the largest harmonic: "+num2str(Va1rms)+"V")
    if (L ~= 0 && C ~= 0 && R ~= 0)
        I_peak = sqrt(sum(I_m.^2)); 
        disp("Peak Load Current: "+num2str(I_peak));
        I_h = sqrt(I_peak^2-I_m(1)^2)/sqrt(2);
        disp("Harmonic load current: "+num2str(I_h));
        THD = sqrt(I_peak^2-I_m(1)^2)/I_m(1);
        disp("THD: "+num2str(THD*100)+"%");
        I_o1 = I_m(1)/sqrt(2);
        disp("Fundamental harmonic current RMS: "+num2str(I_o1));
        Po = I_o1^2*R;
        disp("Fundamental Output Power: "+num2str(Po));
        Is = Po/Vrms;
        disp("Source current: "+num2str(Is));
        Ip = I_peak;
        disp("Peak transistor current: "+num2str(Ip));
        Iqmax = Ip/2;
        disp("Permissible transistor current: "+num2str(Iqmax));
    end
end

