function A=spiral_inductance_analysis(filename,offset)
%==========================================================================
%A function used to extract the inductance and the parasitic capacitance of
%a spiral inductor analysed with Sonnet in a CSV file filename
%==========================================================================
close all
tic
subindex = @(A,r,c) A(r,c); 
M=csvread(filename,offset,0);
frequency_temp=M(:,1);
index_cut=[];
fmax=100;
fmin=0.1;

for i=1:size(M,1)
    if(frequency_temp(i)<fmax)&(frequency_temp(i)>fmin)
        index_cut=[index_cut,i];
    end
end
M=subindex(M,index_cut,:);
frequency=(10^9).*M(:,1);
Re_Y11=M(:,2);
Im_Y11=M(:,3);
Re_Y12=M(:,4);
Im_Y12=M(:,5);
Re_Y21=M(:,6);
Im_Y21=M(:,7);
Re_Y22=M(:,8);
Im_Y22=M(:,9);
Y12=Re_Y12+1i*Im_Y12;
Y11=Re_Y11 +1i*Im_Y11;
Y21=Re_Y21+1i*Im_Y21;
Y22=Re_Y22 +1i*Im_Y22;
L=(1./(2*pi.*frequency)).*imag(1./Y11);
L_diff=(4./(2*pi.*frequency)).*imag(1./(Y11+Y22-Y21-Y12));
[maxi,I_maxi]=max(abs(L));
[maxi_diff,I_maxi_diff]=max(abs(L_diff));
omega_0=frequency(I_maxi)
omega_0_diff=frequency(I_maxi_diff)
L_i=L(1);
L_i_diff=L_diff(1);
C_para=1/(L_i*4*pi*pi*omega_0*omega_0)
C_para_diff=1/(L_i_diff*4*pi*pi*omega_0_diff*omega_0_diff)
F=@(X,xdata)X(1)./(1-(xdata.^2)/(X(2)^2));
% X0=[L_i, omega_0];
% X=lsqcurvefit(F,X0,frequency,L)

figure(1)
subplot(2,1,1)
plot(frequency,L)
title({['L =' num2str((10^9)*L_i) 'nH C =' num2str((10^15)*C_para) 'fF'];['omega_0 = ' num2str(omega_0*10^(-9)) 'Ghz']})
xlabel('frequency (Ghz)')
ylabel('L(omega) (H)')
subplot(2,1,2)
plot(frequency,L_diff)
title({['L_diff =' num2str((10^9)*L_i_diff) 'nH C_diff =' num2str((10^15)*C_para_diff) 'fF'];['omega_0_diff = ' num2str(omega_0_diff*10^(-9)) 'Ghz']})
xlabel('frequency (Ghz)')
ylabel('L_diff(omega) (H)')
% % subplot(2,2,3)
% % plot(frequency,L_diff_bis)
% % title({['L_diff =' num2str((10^9)*L_i_diff_bis) 'nH C_diff =' num2str((10^15)*C_para_diff_bis) 'fF'];['omega_0_diff = ' num2str(omega_0_diff_bis*10^(-9)) 'Ghz']})
% % xlabel('frequency (Ghz)')
% % ylabel('L_diff_bis(omega) (H)')
% % figure(2)
% % hold on
% % plot(frequency,L)
% % frequency_bis=linspace(frequency(1),frequency(size(frequency,1)),1000);
% % plot(frequency_bis,L_i./(1-(frequency_bis.^2)/(omega_0^2)),'-r')
% % plot(frequency_bis,L_i_diff./(1-(frequency_bis.^2)/(omega_0_diff^2)),'-g')
% % hold off
toc
