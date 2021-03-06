function A=treating_sonnet_result(filename,offset)
%==========================================================================
%A function used to extract the resonance frequency, and quality factor (of
%coupling and intrinsec) of a resonator analysed with Sonnet in a CSV file
%filename. The precision should be at least 1.0e-4/1.0e-5 Ghz around the
%resonance to have a good analysis. The resonators should have two ports,
%otherwise one need to change the way to extract the S21 parameter
%==========================================================================
tic
close all
subindex = @(A,r,c) A(r,c); 



M=csvread(filename,offset,0);
frequency_temp=M(:,1).*10^(-9);
% frequency_temp=M(:,1);
u=size(frequency_temp)
index_cut=[];
fmax=20;
fmin=1;
for i=1:size(M,1)
    if(frequency_temp(i)<fmax)&(frequency_temp(i)>fmin)
        index_cut=[index_cut,i];
    end
end
M=subindex(M,index_cut,:);
% [mini,I]=min(M(:,6).^2+M(:,7).^2)
[mini,I]=min(real(M(:,2)).^2+imag(M(:,2)).^2)
frequency_temp=M(:,1).*10^(-9);
% frequency_temp=M(:,1);
omega_0_guess=frequency_temp(I)
bandwidth=0.01;
L=[];
for i=1:size(M,1)
    if(abs(frequency_temp(i)-omega_0_guess)<bandwidth)
        L=[L,i];
    end
end
M=subindex(M,L,:);
v=size(M)
% [min_value,I]=min(M(:,6).^2+M(:,7).^2);
[min_value,I]=min(real(M(:,2)).^2+imag(M(:,2)).^2);
n=size(M,1);
% frequency=M(:,1);
frequency=M(:,1).*10^(-9);
format long
% Re_S11=M(:,2);
% Im_S11=M(:,3);
% Re_S12=M(:,4);
% Im_S12=M(:,5);
% Re_S21=M(:,6);
% Im_S21=M(:,7);
% Re_S22=M(:,8);
% Im_S22=M(:,9);
Re_S21=real(M(:,2));
Im_S21=imag(M(:,2));

[Mbis,Ibis]=max(Re_S21);
[M3,I3]=max(Im_S21)
[M4,I4]=min(Im_S21)

real_bias=0.5*(Re_S21(1)+Re_S21(n))
imag_bias=0.5*(Im_S21(1)+Im_S21(n))
Re_S21=Re_S21./real_bias;
Im_S21=Im_S21/real_bias;

u=1



Q_0_guess=(abs(frequency(I3)-frequency(I4))/omega_0_guess)^(-1)
Q_c_guess=Q_0_guess*((u-Re_S21(I))^(-1))
Q_i_guess=Q_c_guess*Q_0_guess/(abs(Q_c_guess-Q_0_guess))

real_bias-Re_S21(I)
a=Re_S21(1)
c=Re_S21(n)
b=real_bias


% delta_omega_guess=(-0.5)*(Q_0_guess/Q_c_guess)*(frequency(I3)+frequency(I4)-2*omega_0_guess)
delta_omega_guess=(-0.5)*omega_0_guess/Q_0_guess

offset_RE=0.5*(Re_S21(1)+Re_S21(n)-2);
offset_IM=0.5*(Im_S21(1)+Im_S21(n));

% % % F_real=@(X,xdata)offset_RE+(1-(Q_0_guess/Q_c_guess)+4*((Q_0_guess/omega_0_guess)^2).*(xdata-omega_0_guess).*(xdata-omega_0_guess+X))./(1+4*((Q_0_guess/omega_0_guess)^2).*((xdata-omega_0_guess).^2));
% % % F_imag=@(X,xdata)offset_IM+(2*Q_0_guess/omega_0_guess)*((Q_0_guess/Q_c_guess).*(xdata-omega_0_guess)+X)./(1+4*((Q_0_guess/omega_0_guess)^2).*((xdata-omega_0_guess).^2));
% % % X0=[delta_omega_guess];
Qtot=@(X)(X(1)*X(2))/((X(1)+X(2)));

F_tot=@(X,xdata)offset_RE+offset_IM*1i+u-(Qtot([X(1),X(2)])/X(2)-2i*Qtot([X(1),X(2)])*X(3)/omega_0_guess)./(1+2i*Qtot([X(1),X(2)])/omega_0_guess.*(xdata-omega_0_guess))

F_real=@(X,xdata)X(4)+u*(1-(Qtot([X(1),X(2)])/X(2))+4*((Qtot([X(1),X(2)])/omega_0_guess)^2).*(xdata-omega_0_guess).*(u.*(xdata-omega_0_guess)+X(3)))./(1+4*((Qtot([X(1),X(2)])/omega_0_guess)^2).*((xdata-omega_0_guess).^2));
F_imag=@(X,xdata)X(5)+(2*Qtot([X(1),X(2)])/omega_0_guess)*(u*(Qtot([X(1),X(2)])/X(2)).*(xdata-omega_0_guess)+X(3))./(1+4*((Qtot([X(1),X(2)])/omega_0_guess)^2).*((xdata-omega_0_guess).^2));
X0=[Q_i_guess,Q_c_guess,delta_omega_guess,offset_RE,offset_IM];

X=lsqcurvefit(F_imag,X0,frequency,Im_S21)
X=lsqcurvefit(F_real,X,frequency,Re_S21)


omega_0_guess
Q_c=X(2)
Q_c_guess
Q_i=X(1)
Q_i_guess
Q_0=Qtot([X(1), X(2)])
Q_0_guess
% Q_i=Q_c*Q_0/(abs(Q_c-Q_0))
delta_omega=X(3)
delta_omega_guess

figure(2)
subplot(2,2,1)
hold on
plot(F_real(X0,frequency),F_imag(X0,frequency),'-r')
plot(Re_S21,Im_S21)
title(['Q_0 = ', num2str(Q_0_guess), ', Q_i = ', num2str(Q_i_guess)])
xlabel('Re(S21)')
ylabel('Im(S21)')
hold off
subplot(2,2,2)
hold on
plot(frequency,Re_S21)
plot(frequency,F_real(X0,frequency),'-r')
title(['Q_c = ', num2str(Q_c_guess)])
xlabel('{frequency}{(Ghz)}')
ylabel('Re(S21)')
hold off
subplot(2,2,3)
hold on
plot(frequency,Im_S21)
plot(frequency,F_imag(X0,frequency),'-r')
title([' omega_0 = ', num2str(omega_0_guess), ' Ghz'])
xlabel('{frequency}{(Ghz)}')
ylabel('Im(S21)')
hold off
subplot(2,2,4)
hold on
plot(frequency,Im_S21.^2+Re_S21.^2)
plot(frequency,F_real(X0,frequency).^2+F_imag(X0,frequency).^2,'-r')
xlabel('{frequency}{(Ghz)}')
ylabel('Magnitude of S21')
hold off


figure(1)
subplot(2,2,1)
hold on
plot(F_real(X,frequency),F_imag(X,frequency),'-r')
plot(Re_S21,Im_S21)
title(['Q_0 = ', num2str(Q_0), ', Q_i = ', num2str(Q_i)])
xlabel('Re(S21)')
ylabel('Im(S21)')
hold off
subplot(2,2,2)
hold on
plot(frequency,Re_S21)
plot(frequency,F_real(X,frequency),'-r')
title(['Q_c = ', num2str(Q_c)])
xlabel('{frequency}{(Ghz)}')
ylabel('Re(S21)')
hold off
subplot(2,2,3)
hold on
plot(frequency,Im_S21)
plot(frequency,F_imag(X,frequency),'-r')
title([' omega_0 = ', num2str(omega_0_guess), ' Ghz'])
xlabel('{frequency}{(Ghz)}')
ylabel('Im(S21)')
hold off
subplot(2,2,4)
hold on
plot(frequency,Im_S21.^2+Re_S21.^2)
plot(frequency,F_real(X,frequency).^2+F_imag(X,frequency).^2,'-r')
xlabel('{frequency}{(Ghz)}')
ylabel('Magnitude of S21')
hold off

figure(3)
subplot(2,2,1)
hold on
plot(real(F_tot(X0,frequency)),imag(F_tot(X0,frequency)),'-r')
plot(Re_S21,Im_S21)
title(['Q_0 = ', num2str(Q_0_guess), ', Q_i = ', num2str(Q_i_guess)])
xlabel('Re(S21)')
ylabel('Im(S21)')
hold off
subplot(2,2,2)
hold on
plot(frequency,Re_S21)
plot(frequency,real(F_tot(X0,frequency)),'-r')
title(['Q_c = ', num2str(Q_c_guess)])
xlabel('{frequency}{(Ghz)}')
ylabel('Re(S21)')
hold off
subplot(2,2,3)
hold on
plot(frequency,Im_S21)
plot(frequency,imag(F_tot(X0,frequency)),'-r')
title([' omega_0 = ', num2str(omega_0_guess), ' Ghz'])
xlabel('{frequency}{(Ghz)}')
ylabel('Im(S21)')
hold off
subplot(2,2,4)
hold on
plot(frequency,Im_S21.^2+Re_S21.^2)
plot(frequency,real(F_tot(X0,frequency)).^2+imag(F_tot(X0,frequency)).^2,'-r')
xlabel('{frequency}{(Ghz)}')
ylabel('Magnitude of S21')
hold off


toc