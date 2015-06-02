function A=spiral_inductance_modeling(n,Do,s,w)
geometrical_factor=@(S,W)ellipke(sqrt(1-(1/(1+2*(S/W)))^2))/ellipke(1/(1+2*(S/W)))
L=Do;
C=0;
epsilon0=8.85*10^(-12);
epsilonr=1;
epsilon=epsilon0*epsilonr;
for i=1:(n-1)
    C=C+(L-2*(s+w))*10^(-6)*epsilon*geometrical_factor(s,w);
    L=L-2*(s+w);
end
C=C*10^15
