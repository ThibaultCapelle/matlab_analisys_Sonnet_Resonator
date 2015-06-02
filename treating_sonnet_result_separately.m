function A=treating_sonnet_result_separately(filename1,filename2,offset1,offset2)
tic
close all
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
    S_final(:,:,j)=subindex(L(:,:,j),[1,2],[1,2])+D(:,:,j)*subindex(L(:,:,j),[3,4],[1,2])*subindex(L(:,:,j),[3,4],[3,4]);
end   
ReS21=squeeze(real(subindex3(S_final,2,1,:)));
ImS21=squeeze(imag(subindex3(S_final,2,1,:)));
figure(1)
plot(frequency_i,ReS21)
figure(2)
plot(frequency_i,ImS21)

toc