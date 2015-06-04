function A=capacitance_membrane(filename,offset)
%==========================================================================
%Function used to analyse the capacitance measured with Sonnet and saved in
%the CSV file filename and to try to fit the model to the graph
%==========================================================================
close all
tic
hmax=100;%the maximum altitude value used to troncate the data
M=csvread(filename,offset,0);
subindex = @(A,r,c) A(r,c); 
subind=[];
for i=1:size(M,1)
    if(M(i,1)<hmax)
        subind=[subind,i];
    end
end
M=subindex(M,subind,:);%the data is troncated before any analysis
thickness=M(:,1);
capacitance=M(:,2);
n_c=size(capacitance,1)
C_0=capacitance(n_c);%the capacitance without membrane is assumed to be the capacitance with the highest altitude found in the data analysed
epsilon_membrane=9.61;
r=(1-epsilon_membrane)/(1+epsilon_membrane);
unity=ones(n_c,1);
a=2;
s=1;
Nmax=1000;%the highest order in the Fourier decomposition to be used as a model
for i=999:Nmax
    %=========================================================
    %calculation of the model at the order i
    %=========================================================
    N=i
    n=[1:N];
    alpha=(3.141592/a).*(2.*n-1);
    B=((4/a)./alpha).*besselj(0,alpha.*s/2);
    norm=sum(B.*cos(alpha*s/2));
    Capatot=@(X,xdata)(C_0*norm)./sum(kron(unity,B).*cos(kron(unity,alpha)*s/2).*(1+r*exp(-X*2*kron(xdata,alpha))),2);
%     X0=[0.5];
%     X=lsqcurvefit(Capatot,X0,thickness,capacitance)
    figure(1)
    hold on;
    plot(thickness,capacitance,'LineWidth',2)
    % plot(thickness,Capatot(r,thickness),'-r')
    plot(thickness,Capatot(1,thickness))
end
r
C_0
toc