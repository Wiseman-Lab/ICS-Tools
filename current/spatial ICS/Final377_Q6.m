clear all
close all

w1 = 530; %nm
w2 = 531; %nm
v1 = 1/w1; %nm-1
v2 = 1/w2; %nm-1
FWHM1 = 1.1;
FWHM2 = 1.1;
Length = 100000;

while FWHM1 > 1 || FWHM2 > 1
    s = 100; %Sampled every s nm
    Length = Length + 10000; %Sampling length L (nm)
    L = Length/s; %# of data points
    
    b = [0:L-1]*s;
    
    A1 = 1;
    A2 = A1; %Amplitudes are equal
    
    p1 = 0;
    p2 = 0;
    
    y = A1*cos(2*pi*v1*b) + A2*cos(2*pi*v2*b);
    figure(1)
    plot(b,y)
    xlabel('x(nm)')
    ylabel('Amplitude')
    
    %%
    
    Y = fft(y);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = 1/s*(0:(L/2))/L;
    
    figure(2)
    plot(f(:,1:(size(f,2)/100)),P1(:,1:(size(P1,2)/100)))
    xlabel('Wavenumber (nm^-^1)')
    ylabel('Power')
    
    
    W = 1./f;
    figure(3)
    plot(W,P1)
    xlabel('Wavelength (nm)')
    ylabel('Power')
    axis([520 540 0 2])
    
    %% Determining FWHM
    pks = findpeaks(P1,'MinPeakHeight',0.5);
    
    if size(pks,2) == 2
        ind1 = find(P1 == pks(1))
        ind2 = find(P1 == pks(2))
        loc1 = W(1,ind1);
        loc2 = W(1,ind2);
        
        FWHM1 = abs(W(1,ind1+1)-W(1,ind1-1))
        FWHM2 = abs(W(1,ind2+1)-W(1,ind2-1))
    else
        FWHM1 = 2;
        FWHM2 = 2;
    end
end

Length
s