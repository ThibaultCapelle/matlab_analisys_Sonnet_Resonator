close all
Qc=50000
Qi=10
Q=1.0/(1.0/Qc+1.0/Qi)
omega0=10
omega=linspace(0,20,10000)
deltaomega=omega-omega0
S=1-(Q/Qc)./(1+2i*(Q/omega0).*deltaomega)
figure(1)
plot(omega,abs(S))