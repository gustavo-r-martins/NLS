function [L,beta2,beta3,gamma,N,L_D,L_NL] = param(L,beta2,beta3,gamma,P0,T0)
    arguments
        L = 1
        beta2 = 1
        beta3 = 0
        gamma = 0
        P0 = 1
        T0 = 1
    end    
    if (beta3 ~= 0)
        L_D = T0^3/abs(beta3);
        N = sqrt(gamma*P0*T0^3/abs(beta3));                     % soliton order
        beta3 = sign(beta3)/L_D;
    end
    if (beta2 ~= 0)
        L_D = T0^2/abs(beta2);
        N  = sqrt(gamma*P0*T0^2/abs(beta2));                    % soliton order
        beta2 = sign(beta2)/L_D;
    end

    if gamma ~= 0
        L_NL = 1/(gamma*P0);
        gamma = 1/L_NL;
    end

%     if alpha ~= 0
%         Leff = (1 - exp(-alpha*L))/alpha;
%     else
%         Leff = L;
%     end
end