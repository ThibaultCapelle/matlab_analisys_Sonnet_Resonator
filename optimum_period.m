close all
epsilon_INP=9.61
epsilon_SiO2=3.9
epsilon_Si=11.68
h0=0.26e-6
h=0.05e-6
h_resine=2e-6

amin=1e-9
amax=10e-6
a=linspace(amin,amax,(amax-amin)/amin+1);
alpha=pi./a;



reff=-(epsilon_INP^2-1).*(exp(2*pi*h0./a)-1)./(exp(2*pi*h0./a).*(1+epsilon_INP)^2-(epsilon_INP-1)^2);
r1=(1-epsilon_SiO2)/(1+epsilon_SiO2)
t1=2/(1+epsilon_SiO2)
r1b=(epsilon_SiO2-1)/(1+epsilon_SiO2)
t1b=2*epsilon_SiO2/(1+epsilon_SiO2)
r2=(epsilon_SiO2-epsilon_Si)/(epsilon_Si+epsilon_SiO2)
rsub=r1.*ones(1,size(a,2))+t1*t1b*r2.*exp(-2*pi*h_resine./a)./(1-r1b*r2.*exp(-2*pi*h_resine./a));






factor=(1-rsub.*reff.*exp(-2*alpha.*h))./((1+rsub).*(1+reff.*exp(-2*alpha.*h)));
G=abs(pi*reff.*(1+rsub).*exp(-2*pi*h./a)./(a.*(1-reff.*rsub.*exp(-2*pi*h./a)).*(1+reff.*exp(-2*pi*h./a))))

a=a(not(isnan(G)))
G=G(not(isnan(G)))

size(G,2)
G_max=G(1)
I_max=1
for i=1:size(G,2)
    if G(i)>=G_max
        I_max=i;
        G_max=G(i)
    end
end
a_ideal=a(I_max)
        

% u=max(G,[],'omitnan')


figure(1)
plot(a,G)
xlabel('a(m)')
ylabel('Coupling (u.a.)')
title(['maximum coupling rate for a = ',num2str(a_ideal/(1e-9)),' nm'])
h=legend(strvcat('membrane dielectric constant : 3.6','membrane thickness : 0.26 microns','height of membrane : 50 nm','thickness of oxyde : 2 microns'));
h.FontSize=12