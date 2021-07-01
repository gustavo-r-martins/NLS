% This code solves the NLS equation with the split-step method
%     idu/dz - sgn(beta2)/2 d^2u/d(tau)^2 + N^2*|u|^2*u = 0
%     idu/dz = beta2/2 d^2u/d(tau)^2 + i beta3/6 d^3u/d(tau)^3
% Written by Govind P. Agrawal in March 2005 for the NLFO book

% Specify input parameters
clear;

L = 10.6;                 % Fiber length (in units of L_D);
beta2 = -5.23e-27;      % dispersion: 1 for normal, -1 for anomalous;
beta3 = 4.27e-41;       % TOD parameter: 1 for right, -1 for left;
gamma = 18.4e-3;        % nonlinear parameter 
alpha = 0;              % fibre loss (~0.2 dB/km)
S = 0;                  % self-stepening
tauR = 0;               % Raman scatering

mshape = 0;             % m = 0 for sech, m > 0 for super-Gaussian;
chirp0 = 0;             % pulse chirp (default value C = 0)

P0 = 26.3;                 % input pulse peak power
TFWHM = 1.1e-12;
lambda0 = 1550e-6;           % central wavelength
c = physconst('LightSpeed');

if mshape==0
    T0 = TFWHM/(2*log(1+sqrt(2)));
else
    T0 = TFWHM/(2*sqrt(log(2)));
end
% T0 = sqrt(N^2*abs(beta2)/gamma/P0); % input pulse width

[L,beta2,beta3,gamma,N,L_D,L_NL] = param(L,beta2,beta3,gamma,P0,T0);

%---set simulation parameters
nt = 2048; Tmax = 32;                                                      % FFT points ad window size
step_num = 10e3;                                                           % No. of steps
deltaz = L/step_num;                                                       % step size in z
dtau = (2*Tmax)/nt;                                                        % step size in tau
fs = 1/dtau;

%---tau and omega arrays
tau = (-nt/2:nt/2-1)*dtau;                                                 % temporal grid
omega = (2*pi*fs/nt)*[(0:nt/2-1) (-nt/2:-1)];                              % frequency grid

%---Input Field profile
if mshape==0
    uu = sech(tau).*exp(-0.5i*chirp0*(tau).^2);                              %soliton
else
    uu = sqrt(P0)*exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape));                         % super-Gaussian
end

%---Plot input pulse shape and spectrum
temp = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi);                           %spectrum
figure(1); hold on;
subplot(2,1,1);
    plot(tau*T0, P0*abs(uu).^2, '--k');
    hold on;
    axis([-2e-12 2e-12 0 Inf]);
    xlabel('Time');
    ylabel('Potência (W)');
%     legend('Input')
    title('Evolução do Pulso Simulado com NLSE');
subplot(2,1,2);
    plot(fftshift(omega)/(2*pi)/T0, db(P0*abs(temp).^2), '--k');
    hold on;
    axis([-Inf Inf -40 40]);
    xlabel('(w-w0)T_0');
    ylabel('Potência Espectral (dB)');

% Fiber Hamiltonians
% (w-w0) -> i(d/dt)
%---store dispersive phase shifts to speedup code
D = 1i*beta2*omega.^2/factorial(2);
D = 1i*beta3*omega.^3/factorial(3) + D;

D = exp(D*deltaz);                                                          % phase factor 
NL = 1i*gamma*deltaz;
% NL = 1i*N^2*deltaz;                                                       % nonlinear phase factor


% scheme: 1/2N -> D -> 1/2N; first half step nonlinear
% Cavity
uu = ring(uu,D,step_num,NL,S,omega,tauR);


temp = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi);                            % Final spectrum

%---Plot output pulse shape and spectrum
subplot(2,1,1)
    plot(tau*T0, P0*abs(uu).^2,'k')
    legend('Input','Output')
subplot(2,1,2)
    plot(fftshift(omega)/(2*pi)/T0, db(P0*abs(temp).^2),'k')





