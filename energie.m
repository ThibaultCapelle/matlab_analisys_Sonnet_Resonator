close all;
N=100;
u=1:N;
M=1000;
L=linspace(0,5,M);
B=@(lambda,n)(1./(2*n-1)).*besselj(0,pi*lambda.*(2*n-1)./2);
E=@(lambda,n)(B(lambda,n).^2)*(2*n-1);
energy=zeros(N,M);
for i=1:N
    energy(i,:)=E(L,i);
end
figure(1)
hold on
for i=1:N
    plot(L,energy(i,:))
end
hold off
figure(2)
plot(L,energy(1,:)./sum(energy(2:N,:),1))
xlabel('s/a')
ylabel('predominance of the first component')