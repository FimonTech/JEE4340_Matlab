function construct_signal(x,L,Fs,t,N,name)
    X = fft(x);         % this is where DSP magic happens
    P2 = abs(X/L);  
    P1 = P2(1:floor(L/2+1));
    P1(2:end-1) = 2*P1(2:end-1); % removing symmetrical frequencies at the end
    a0 = real(P1(1));  % a0 = 2 * mean(x)
    an = zeros(1,floor(L/2+1));  % Depends on how many harmonics we'll have here
    bn = zeros(1,floor(L/2+1));  % Ill just keep the array big enough, it's fine
    
    f = Fs*(0:floor(L/2))/L;
    figure()
    [psor,lsor] = findpeaks(P1,f,'SortStr','descend');
    plot(f,P1) 
    title("Single-Sided Amplitude Spectrum of X(t)")
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    text(lsor+.02,psor,num2str((1:numel(psor))'))
    
    v2 = a0 * ones(size(t));
    v2_wo1 = a0 * ones(size(t)); % Waveform without the first harmonic
    %{
    disp("Harmonic number: 0"); % this is to check which harmonic worked
    disp("Harmonic ampltiude: "+num2str(a0)); % this is to check which harmonic worked
    disp("Harmonic frequency: 0")
    %}
    i = 2; % Calculating harmonics through this crap
    n = 1; %iterator
    harmonics = ones(1,10);
    while n <= N % Count to 10 harmonics
        an(i) = 2*real(X(i)/L);
        bn(i) = -2*imag(X(i)/L);
        if (P1(i) > 0.06) % tolerance
            v2 = v2 + an(i).*cos(2*pi*f(i).*t)+bn(i).*sin(2*pi*f(i).*t); % we only care about non-zero harmonics
            harmonics(n) = sqrt(an(i)^2+bn(i)^2); % Get the sum of harmonics (needed for RMS calculation)
            if (n > 1)
               v2_wo1 = v2_wo1 + an(i).*cos(2*pi*f(i).*t)+bn(i).*sin(2*pi*f(i).*t);
            end
            n = n + 1;
            %{
            disp(an(i));
            disp(bn(i));
            disp("Harmonic number: "+num2str(i)); % this is to check which harmonic worked
            disp("Harmonic amplitude: "+num2str(P1(i))); % this is to check which harmonic worked
            disp("Harmonic frequency: "+num2str(f(i)))
            %}
        end
        i = i + 1;
    end
    plot(t*1000,v2);
    hold on 
    plot(t*1000,v2_wo1,'--',"Color","#EDB120")
    hold off
    title("Reconstructed signal of "+name)
    xlabel("Time, ms")
    ylabel("V(t), V")
    legend(["Reconstructed V_0(t)","Reconstructed V_0(t) w/o 1st harmonic"])
    grid on
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
end

