function [Mu U]=respspec(a,dt,f,z,forcing)

% D.Melgar 02/2012
%
% Compute displacement response spectra for acceleration or displacement
% forcing by solving the single dgree of freedom oscillator ODE using 
% second order accurate finite differences
%
% IN
% a ~ Time series, may be acceleration or displacement time series
% dt ~ Sampling interval in seconds
% f ~ Vector of frequencies at which responses will be computed, usually a 
%     vector with some 200-500 elements works nicely but more is possible 
%     if you are willing to wait for the computations.
% z ~ Damping in decimals i.e. 0.01=1%
% forcing ~ 'a' if forcing time series is acceleration and 'd' if it is
%           displacement. I did not include velocity because I am lazy but shoot me an 
%			email and we'll talk about it
%
% OUT
% Mu ~ Displacement response spectra
% U ~ Displacement time histories for the oscillator at all input frequencies
%
% CITING THIS CODE
%
% Please reference
% D.Melgar, Bock, Y., Snachez, D. & Crowell, B.W. (2013). On Robust and Automated Baseline
% Corrections for Strong Motion Seismology, J. Geophys. Res, 118(3), 1177-1187, DOI: 10.1002
% /jgrb.50135.

%Make row vector
N=size(a);
if N(1)>N(2)
    a=a';
end
N=max(N);
%Zero pad for long period spectra
a=[zeros(1,N) a zeros(1,N)];
i=N+1:N*2;
%How many frequencies?
Nf=max(size(f));
%Set up the run
w=2*pi*f;%omega. Here w has Nf values
U=zeros(Nf,N);
u=zeros(1,N*3);
%Compute spectra
for k=1:Nf;
    wc=w(k);%omega current, wc. This will retrieve omega(w)values one at a time.
    %Solve
    c=(1/(dt^2)+(z*wc)/dt);%c is just a bunch of constants
    if forcing=='a'
        for k2=3:N*3
            %Accleration forcing
            u(k2)=(1/c)*((2*u(k2-1)-u(k2-2))/(dt^2)+((z*wc*u(k2-2))/dt)-u(k2-1)*wc^2-a(k2-1));
        end
    elseif forcing=='d'
        for k2=3:N*3
            %Displacement forcing
            u(k2)=(1/c)*(((2*u(k2-1)-u(k2-2))/(dt^2)+((z*wc*u(k2-2))/dt)-u(k2-1)*wc^2+(-a(k2)+2*a(k2-1)-a(k2-2))/(dt^2)));
        end
    end
    %Save for output
    U(k,:)=u(i);
    Mu(k)=max(abs(u(i)));
end


