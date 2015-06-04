parameter=[(1/2) 1];
geometric_factor=ellipke(sqrt(1-1./(1+1.*(parameter)).^2))./ellipke(1./(1+1.*(parameter)));
geometric_approx=1./log(1+(1./parameter)+sqrt((1+1./parameter).^2-1));
% for i=1:size(parameter)
%     k=1/(1+2*parameter(i));
%     k2=sqrt(1-k^2);
%     if(k<=(1/sqrt(2)))
%         geometric_approx(i)=log(2*(1+sqrt(k2))/(1-sqrt(k2)))/pi;
%     else
%         geometric_approx(i)=pi/log(2*(1+sqrt(k))/(1-sqrt(k)));
%     end
% end
% figure(1)
% hold on
% plot(parameter,geometric_factor)
% % plot(parameter,geometric_approx,'-r')
% title('geometric factor value as a function of w/s') 
u=geometric_factor(2)/geometric_factor(1)
v=geometric_factor(1)
w=geometric_factor(2)
z=geometric_approx(2)/geometric_approx(1)
% v=ellipke(sqrt(1-1./(1+1.*(1/2)).^2))./ellipke(1./(1+1.*(1/2)))
% w=ellipke(sqrt(1-1./(1+2.*(1/3)).^2))./ellipke(1./(1+2.*(1/3)))
% z=w/v