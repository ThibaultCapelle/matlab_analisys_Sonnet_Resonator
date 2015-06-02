function A=treating_sonnet_result_separately_compare(filename1,filename2,filename3,offset1,offset2,offset3)
tic
close all
subindex1 = @(A,r) A(r);
subindex = @(A,r,c) A(r,c); 
subindex3 = @(A,r,c,l) A(r,c,l); 
M=csvread(filename1,offset1,0);
N=csvread(filename2,offset2,0);
m=size(M,1)
n=size(N,1)
k=min(m,n);
M=subindex(M,[1:k],:);
N=subindex(N,[1:k],:);
frequency_c=M(:,1);
frequency_i=N(:,1);

Re_S11_c=M(:,2);
Im_S11_c=M(:,3);
Re_S12_c=M(:,4);
Im_S12_c=M(:,5);
Re_S21_c=M(:,6);
Im_S21_c=M(:,7);
Re_S22_c=M(:,8);
Im_S22_c=M(:,9);

Re_S11_i=N(:,2);
Im_S11_i=N(:,3);
Re_S12_i=N(:,4);
Im_S12_i=N(:,5);
Re_S13_i=N(:,6);
Im_S13_i=N(:,7);
Re_S14_i=N(:,8);
Im_S14_i=N(:,9);
Re_S21_i=N(:,10);
Im_S21_i=N(:,11);
Re_S22_i=N(:,12);
Im_S22_i=N(:,13);
Re_S23_i=N(:,14);
Im_S23_i=N(:,15);
Re_S24_i=N(:,16);
Im_S24_i=N(:,17);
Re_S31_i=N(:,18);
Im_S31_i=N(:,19);
Re_S32_i=N(:,20);
Im_S32_i=N(:,21);
Re_S33_i=N(:,22);
Im_S33_i=N(:,23);
Re_S34_i=N(:,24);
Im_S34_i=N(:,25);
Re_S41_i=N(:,26);
Im_S41_i=N(:,27);
Re_S42_i=N(:,28);
Im_S42_i=N(:,29);
Re_S43_i=N(:,30);
Im_S43_i=N(:,31);
Re_S44_i=N(:,32);
Im_S44_i=N(:,33);


S11_c=Re_S11_c +1i*Im_S11_c;
S12_c=Re_S12_c +1i*Im_S12_c;
S21_c=Re_S21_c +1i*Im_S21_c;
S22_c=Re_S22_c +1i*Im_S22_c;

S11_i=Re_S11_i +1i*Im_S11_i;
S12_i=Re_S12_i +1i*Im_S12_i;
S21_i=Re_S21_i +1i*Im_S21_i;
S22_i=Re_S22_i +1i*Im_S22_i;
S13_i=Re_S13_i +1i*Im_S13_i;
S14_i=Re_S14_i +1i*Im_S14_i;
S23_i=Re_S23_i +1i*Im_S23_i;
S24_i=Re_S24_i +1i*Im_S24_i;
S33_i=Re_S33_i +1i*Im_S33_i;
S34_i=Re_S34_i +1i*Im_S34_i;
S43_i=Re_S43_i +1i*Im_S43_i;
S44_i=Re_S44_i +1i*Im_S44_i;
S31_i=Re_S31_i +1i*Im_S31_i;
S32_i=Re_S32_i +1i*Im_S32_i;
S41_i=Re_S41_i +1i*Im_S41_i;
S42_i=Re_S42_i +1i*Im_S42_i;
L=zeros(4,4,k);
S_final=zeros(2,2,k);
D=zeros(2,2,k);
u=8
for j=1:k
    L(:,:,j)=[S11_i(j) S12_i(j) S13_i(j) S14_i(j); S21_i(j) S22_i(j) S23_i(j) S24_i(j) ; S31_i(j) S32_i(j) S33_i(j) S34_i(j) ; S41_i(j) S42_i(j) S43_i(j) S44_i(j)];
    D(:,:,j)=([S11_c(j) S12_c(j) ; S21_c(j) S22_c(j)]^(-1)-[S33_i(j) S34_i(j) ; S43_i(j) S44_i(j)])^(-1);
    S_final(:,:,j)=subindex(L(:,:,j),[1,2],[1,2])+subindex(L(:,:,j),[3,4],[3,4])*D(:,:,j)*subindex(L(:,:,j),[3,4],[1,2]);
end   
ReS21=squeeze(real(subindex3(S_final,2,1,:)));
ImS21=squeeze(imag(subindex3(S_final,2,1,:)));


M=csvread(filename3,offset3,0);
frequency_temp=M(:,1);
[mini,I]=min(M(:,6).^2+M(:,7).^2);
omega_0_guess=M(I,1);
bandwidth=0.5;

u=[size(M,1) size(frequency_i,1) size(frequency_c,1)]
L=[];
L_i=[];
L_c=[];
for i=1:size(frequency_i,1)
    if(abs(frequency_i(i)-omega_0_guess)<bandwidth)
        L_i=[L_i,i];
    end
end

ReS21=subindex1(ReS21,L_i);
ImS21=subindex1(ImS21,L_i);
for i=1:size(frequency_c,1)
    if(abs(frequency_c(i)-omega_0_guess)<bandwidth)
        L_c=[L_c,i];
    end
end
for i=1:size(M,1)
    if(abs(frequency_temp(i)-omega_0_guess)<bandwidth)
        L=[L,i];
    end
end
frequency_c=subindex1(frequency_c,L_c);
frequency_i=subindex1(frequency_i,L_i);
result=size(frequency_c)
M=subindex(M,L,:);
[min_value,I]=min(M(:,6).^2+M(:,7).^2);
n=size(M,1);
frequency=M(:,1);
format long
Re_S11=M(:,2);
Im_S11=M(:,3);
Re_S12=M(:,4);
Im_S12=M(:,5);
Re_S21=M(:,6);
Im_S21=M(:,7);
Re_S22=M(:,8);
Im_S22=M(:,9);


[Mbis,Ibis]=max(Re_S21);
[M3,I3]=max(Im_S21);
[M4,I4]=min(Im_S21);
Q_0_guess=((frequency(I3)-frequency(I4))/omega_0_guess)^(-1)
Q_c_guess=Q_0_guess*((1-Re_S21(I))^(-1))

delta_omega_guess=(-0.5)*(Q_0_guess/Q_c_guess)*(frequency(I3)+frequency(I4)-2*omega_0_guess)

F_real=@(X,xdata)(Re_S21(1)-1)+(1-(X(1)/X(2))+4*((X(1)/omega_0_guess)^2).*(xdata-omega_0_guess).*(xdata-omega_0_guess+X(3)))./(1+4*((X(1)/omega_0_guess)^2).*((xdata-omega_0_guess).^2));
F_imag=@(X,xdata)Im_S21(1)+(2*X(1)/omega_0_guess)*((X(1)/X(2)).*(xdata-omega_0_guess)+X(3))./(1+4*((X(1)/omega_0_guess)^2).*((xdata-omega_0_guess).^2));
X0=[Q_0_guess,Q_c_guess,delta_omega_guess];

X=lsqcurvefit(F_real,X0,frequency,Re_S21)
omega_0_guess
Q_c=X(2)
Q_0=X(1)
Q_i=Q_c*Q_0/(Q_c-Q_0)
delta_omega=X(3)

figure(1)
subplot(2,2,1)
hold on
plot(F_real(X,frequency),F_imag(X,frequency),'-r')
plot(Re_S21,Im_S21)
hold off
subplot(2,2,2)
hold on
plot(frequency,Re_S21)
plot(frequency,F_real(X,frequency),'-r')
hold off
subplot(2,2,3)
hold on
plot(frequency,Im_S21)
plot(frequency,F_imag(X,frequency),'-r')
hold off


figure(3)
hold on
plot(frequency,Re_S21,'-r')
plot(frequency_i,ReS21)
hold off
figure(2)
hold on
plot(frequency,Im_S21,'-r')
plot(frequency_i,ImS21)
hold off

% % % S21=Re_S21+1i*Im_S21;
% % % 
% % % Y0_guess=0.5*(S21(1)+S21(n))
% % % Y0_real_guess=real(Y0_guess)
% % % Y0_imag_guess=imag(Y0_guess)
% % % Y1_guess=S21(Ibis)-Y0_guess;
% % % Y1_imag_guess=imag(Y1_guess);
% % % Y1_real_guess=real(Y1_guess);
% % % S21_center_guess=S21-Y0_guess;
% % % delta_omega=frequency-omega_0_guess;
% % % index_temp=Ibis;
% % % for i=Ibis:n
% % %     index_temp=index_temp+1;
% % %     if(abs(S21_center_guess(i))<0.5*abs(Y1_guess))
% % %         break
% % %     end
% % % end
% % % bw=delta_omega(index_temp)
% % % Q_i_guess=omega_0_guess/bw;
% % % 
% % % 
% % % arg_Y1_guess=atan2(imag(S21(Ibis)-Y0_guess),real(S21(Ibis)-Y0_guess));
% % % scale_guess=abs(S21(Ibis)-Y0_guess);
% % % 
% % % F_real_ter=@(X,xdata)X(1)+X(3)*(cos(X(4))+sin(X(4)).*((xdata-omega_0_guess)./X(5)))./(1+(((xdata-omega_0_guess)./X(5))).^2);
% % % F_imag_ter=@(X,xdata)X(2)+X(3)*(sin(X(4))-cos(X(4)).*((xdata-omega_0_guess)./X(5)))./(1+(((xdata-omega_0_guess)./X(5))).^2);
% % % 
% % % 
% % % X0_ter=[Y0_real_guess,Y0_imag_guess,scale_guess,arg_Y1_guess,bw]
% % % 
% % % Xter=lsqcurvefit(F_real_ter,X0_ter,frequency,Re_S21);
% % % 
% % % figure(12)
% % % subplot(2,1,1)
% % % hold on
% % % plot(frequency,Re_S21,'-g')
% % % plot(frequency,F_real_ter(Xter,frequency),'-r')
% % % hold off
% % % subplot(2,1,2)
% % % hold on
% % % plot(frequency,Im_S21,'-g')
% % % plot(frequency,F_imag_ter(Xter,frequency),'-r')
% % % hold off
% % % 
% % % figure(13)
% % % hold on
% % % plot(Re_S21,Im_S21,'-g')
% % % plot(F_real_ter(Xter,frequency),F_imag_ter(Xter,frequency),'-r')
% % % hold off
% % % 
% % % figure(11)
% % % hold on
% % % plot(frequency,abs(S21),'-g')
% % % plot(frequency,F_imag_ter(Xter,frequency).^2+F_real_ter(Xter,frequency).^2,'-r')
% % % hold off
% % % 
% % % Q=omega_0_guess/Xter(5)

toc