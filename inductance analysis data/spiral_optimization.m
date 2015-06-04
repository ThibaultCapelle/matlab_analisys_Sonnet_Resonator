%==========================================================================
%a script used to exploit the semi-empirical model of parasitic
%capacitance and inductance of spiral inductors to optimize them and see
%the main trends of variation of L/C with rho, when we ask for a specific D
%and L.
%==========================================================================
close all;
tic
s=1;
w=1;
Do=30;
% L_ideal=20;
Dmin=100;%the minimal diameter used for the sweep
Dmax=100;%the maximal diameter used for the sweep
Lmin=30;%the minimal inductance used for the sweep
Lmax=30;%the maximal inductance used for the sweep
sweep=linspace(Dmin,Dmax,Dmax-Dmin+1);
sweep2=linspace(Lmin,Lmax,Lmax-Lmin+1);
n=size(sweep,2);
m=size(sweep2,2);
result=zeros(n,m);
count=0;
max=0;
rho=0;
turns=0;
D=0;
for i=1:n
    for j=1:m      
        percentage_of_time=((i-1)*(Lmax-Lmin+1)+j)/((Dmax-Dmin+1)*(Lmax-Lmin+1))
        toc
        Do=sweep(i);
        L_ideal=sweep2(j);
        inductance=@(X)(9.9*10^(-3)).*(4*s.*(-3.*((Do.*X)./((1+X).*(s+w)))+1)+4.*(((Do.*X)./((1+X).*(s+w)))-1).*(((Do.*X)./((1+X).*(s+w)))-1).*(s+w)./X+4*(((Do.*X)./((1+X).*(s+w)))-1).*w./X)-12.723;
        capacitance=@(X)0.24*((Do.*X)./((1+X).*(s+w)))*(s+w).*X.^(-0.8);
        rapport=@(X)inductance(X)./capacitance(X);
        x=sym('x','real');
        res=vpasolve(inductance(x)==L_ideal,x);
        res=res(res<=1);
        res=res(imag(res)==0);
        res=res(res>=((s+w)/(Do-(s+w))));
        if(size(res)==1)
            count=count+1;
            rho_ideal=res;
            N=Do*rho_ideal/((s+w)*(1+rho_ideal));
            L=inductance(rho_ideal);
            C=capacitance(rho_ideal);
            ratio=L/C;
            if(ratio>max)
                max=ratio;
                turns=N;
                L_optimal=L_ideal;
                rho=rho_ideal;
                D=sweep(i);
            end
            result(i,j)=ratio;
        
        end
    end
end
max
turns
rho
D
L_optimal
inductance=@(X)(9.9*10^(-3)).*(4*s.*(-3.*((D.*X)./((1+X).*(s+w)))+1)+4.*(((D.*X)./((1+X).*(s+w)))-1).*(((D.*X)./((1+X).*(s+w)))-1).*(s+w)./X+4*(((D.*X)./((1+X).*(s+w)))-1).*w./X)-12.723;
capacitance=@(X)0.24*((D.*X)./((1+X).*(s+w)))*(s+w).*X.^(-0.8);
rapport=@(X)inductance(X)./capacitance(X);
sweep3=linspace(0.01,1,1000);
figure(3)
plot(sweep3,rapport(sweep3))
figure(2)
surf(sweep2,sweep,result)
title('L/C maximum in function of L and Do')
ylabel('Do in microns')
xlabel('L in nanoHenri')
toc
% inductance=@(X)(9.9*10^(-3)).*(4*s.*(-3.*((Do.*X)./((1+X).*(s+w)))+1)+4.*(((Do.*X)./((1+X).*(s+w)))-1).*(((Do.*X)./((1+X).*(s+w)))-1).*(s+w)./X+4*(((Do.*X)./((1+X).*(s+w)))-1).*w./X)-12.723;
% % capacitance=@(X)0.24*X(1)*(s+w)*X(2)^(-0.8);
% capacitance=@(X)0.24*((Do.*X)./((1+X).*(s+w)))*(s+w).*X.^(-0.8);
% % rapport=@(X)-inductance(X)/capacitance(X);
% rapport=@(X)inductance(X)./capacitance(X);
% x=sym('x','real');
% res=vpasolve(inductance(x)==L_ideal,x);
% res=res(res<=1);
% rho_ideal=res(res>=((s+w)/(Do-(s+w))))
% N=Do*rho_ideal/((s+w)*(1+rho_ideal))
% L=inductance(rho_ideal)
% C=capacitance(rho_ideal)
% ratio=L/C
% rho=linspace(0.01,1,1000);
% figure(1)
% subplot(1,2,1)
% plot(rho,rapport(rho),'-r')
% subplot(1,2,2)
% plot(rho,inductance(rho),'-r')
% % A=[-1 0;1 0;0 1;0 -1];
% % N=75;
% % M=100;
% % result1=zeros(N,M);
% % result2=zeros(N,M);
% % for i=1:N
% %     for j=1:M      
% %         b=[-1;i;1;-(j/M)];
% %         [x,fval] = fmincon(rapport,Xtest,A,b);
% %         result1(i,j)=x(1);
% %         result2(i,j)=x(2);
% %     end
% % end
% % figure(1)
% % surf(result1)
% % figure(2)
% % surf(result2)
% % b=[-1;9;1;-0.1];
% % [x,fval] = fmincon(rapport,Xtest,A,b)
% % inductance(x)