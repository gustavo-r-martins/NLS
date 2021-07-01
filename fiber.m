function A = fiber(U,D,step,gamma,S,omega,tauR)
    arguments
       U
       D
       step
       gamma = 0
       S = 0
       omega = 0
       tauR = 0
    end
    % A = ring(U,D,step,NL)
    %
    % U is the input field profile
    % D is the dispersion hamiltonian
    % NL is the nolinear hamiltonian
    % step is the number of z-steps
    % gamma is nonlinear parameter
    % S is self-stepening parameter
    % Tr is Raman scatering term
    % The function uses FFT, thus the vectors must be 2^n to improvements.
    
    NL = abs(U).^2 + ...
         1i*S*conj(U).*ifft(1i*(omega).*fft(U))+...
        (1i*S - tauR)*ifft(1i*(omega).*fft(abs(U).^2));
    
    A = U.*exp(NL.*gamma/2); % note hhz/2
    for n=1:step
        f_temp = ifft(A).*D;
        U = fft(f_temp);
        NL = abs(U).^2 + ... 
             1i*S*conj(U).*ifft(1i*(omega).*fft(U))+...
            (1i*S - tauR)*ifft(1i*(omega).*fft(abs(U).^2));
        A = U.*exp(NL.*gamma);
    end
    NL = abs(U).^2 + ... 
         1i*S*conj(U).*ifft(1i*(omega).*fft(U)) + ...
        (1i*S - tauR)*ifft(1i*(omega).*fft(abs(U).^2));
    U = A.*exp(-NL.*gamma/2); % Final field
    A = U;
end





