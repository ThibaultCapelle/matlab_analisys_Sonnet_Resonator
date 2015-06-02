function A=treating_sonnet_result_capa(filename1,filename2,offset1,offset2)
tic
close all

M=csvread(filename1,offset1,0);
N=csvread(filename2,offset2,0);
m=size(M,1)
n=size(N,1)
frequency_c=M(:,1);
frequency_i=N(:,1);

Re_Z11_c=M(:,2);
Im_Z11_c=M(:,3);
Re_Z12_c=M(:,4);
Im_Z12_c=M(:,5);
Re_Z21_c=M(:,6);
Im_Z21_c=M(:,7);
Re_Z22_c=M(:,8);
Im_Z22_c=M(:,9);

Re_Z11_i=N(:,2);
Im_Z11_i=N(:,3);
Re_Z12_i=N(:,4);
Im_Z12_i=N(:,5);
Re_Z21_i=N(:,10);
Im_Z21_i=N(:,11);
Re_Z22_i=N(:,12);
Im_Z22_i=N(:,13);

Z11_c=Re_Z11_c +1i*Im_Z11_c;
Z12_c=Re_Z12_c +1i*Im_Z12_c;
Z21_c=Re_Z21_c +1i*Im_Z21_c;
Z22_c=Re_Z22_c +1i*Im_Z22_c;

Z11_i=Re_Z11_i +1i*Im_Z11_i;
Z12_i=Re_Z12_i +1i*Im_Z12_i;
Z21_i=Re_Z21_i +1i*Im_Z21_i;
Z22_i=Re_Z22_i +1i*Im_Z22_i;

det_c=Z11_c.*Z22_c-Z12_c.*Z21_c;
det_i=Z11_i.*Z22_i-Z12_i.*Z21_i;

C=(-1*10^(3))*(2*pi.*frequency_c.*imag(det_c./Z21_c)).^(-1);

L=imag(det_i./Z21_i)./(2*pi.*frequency_i);
% [M,I]=min(Re_S21.^2+Im_S21.^2)
% omega_0_guess=frequency(I)
% 
% [Mbis,Ibis]=max(Re_S21)
% [M3,I3]=max(Im_S21)
% [M4,I4]=min(Im_S21)
% Q_0_guess=((frequency(I3)-frequency(I4))/omega_0_guess)^(-1)
% Q_c_guess=Q_0_guess*((1-Re_S21(I))^(-1))
% 
% delta_omega_guess=(-0.5)*(Q_0_guess/Q_c_guess)*(frequency(I3)+frequency(I4)-2*omega_0_guess)
% 
% F_real=@(X,xdata)(Re_S21(1)-1)+(1-(X(1)/X(2))+4*((X(1)/omega_0_guess)^2).*(xdata-omega_0_guess).*(xdata-omega_0_guess+X(3)))./(1+4*((X(1)/omega_0_guess)^2).*((xdata-omega_0_guess).^2))
% F_imag=@(X,xdata)Im_S21(1)+(2*X(1)/omega_0_guess)*((X(1)/X(2)).*(xdata-omega_0_guess)+X(3))./(1+4*((X(1)/omega_0_guess)^2).*((xdata-omega_0_guess).^2))
% X0=[Q_0_guess,Q_c_guess,delta_omega_guess];
% 
% X=lsqcurvefit(F_real,X0,frequency,Re_S21)
% format long
% omega_0_guess
% Q_c=X(2)
% Q_0=X(1)
% Q_i=Q_c*Q_0/(Q_c-Q_0)
[M,I]=min(abs((1.0/(2*pi))*1./sqrt((10^(-21))*L.*C)-frequency_i));
omega_resonnance=frequency_i(I)
figure(1)
subplot(2,1,1)
% hold on
% plot(F_real(X,frequency),F_imag(X,frequency),'-r')
plot(frequency_c,C)
xlabel({'frequency','(Ghz)'})
ylabel({'Capacity','(pF)'})
subplot(2,1,2)
plot(frequency_i,L)
xlabel({'frequency','(Ghz)'})
ylabel({'Inductance','(nH)'})
% hold off
% subplot(2,2,2)
% hold on
% plot(frequency,Re_S21)
% plot(frequency,F_real(X,frequency),'-r')
% hold off
% subplot(2,2,3)
% hold on
% plot(frequency,Im_S21)
% plot(frequency,F_imag(X,frequency),'-r')
% hold off
figure(2)
hold on
plot(frequency_i,(1.0/(2*pi))*1./sqrt((10^(-21))*L.*C))
xlabel({'frequency','(Ghz)'})
ylabel({'1/sqrt(LC)','(Hz)'})
plot(frequency_i,frequency_i*10^9)
hold off
toc