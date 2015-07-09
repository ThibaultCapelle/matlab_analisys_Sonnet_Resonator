
hmini=0.001;
hmaxi=0.1;
hstep=0.0001;
filename='resonance_frequency_optimal.csv';
% filename='result_F1.csv';
N=(hmaxi-hmini)/hstep+1;
h=linspace(hmini,hmaxi,N);

if (exist(filename,'file')>0)
    M=csvread(filename);
    Hread=M(:,1);
    factor=10/min(min(M(:,1)),min(min(h),hstep));
    h_temp=setdiff(int64(factor*h),int64(factor*Hread));
    h=double(h_temp)/factor;
    N=size(h,2);
end

inductor_ori=SonnetProject('design spiral try_membrane.son');
inductor=inductor_ori.clone();
inductor.saveAs('inductorbis.son');

result=zeros(N,2);
tic
start=cputime;
for i=1:N
    percentage=(i-1)/N
    h(i)
    toc
    capa_ori=SonnetProject('capacitance only.son');
    capa=capa_ori.clone();
    index=capa.findVariableIndex('thickness');
    capa.modifyVariableValue('thickness',h(i));
    capa.saveAs('capabis.son');


    netlist=SonnetProject();
    netlist.initializeNetlist();
    netlist.addAbsFrequencySweep(2,8);
    netlist.addProjectFileElement('capabis.son',[4,3],1);
    netlist.addProjectFileElement('inductorbis.son',[3,4,1,2],1);
    netlist.addTouchstoneOutput;
    netlist.saveAs('netlistbis.son');
    netlist.simulate('-t');

    [freq, data, freq_noise, data_noise, Zo] = SXPParse('netlistbis.s2p');

    S21=squeeze(data(2,1,:));
    [~,I]=max(diff(abs(S21)));
    omega_0=freq(I)
    result(i,:)=[h(i),omega_0];
    duration=cputime-start
    estimated_time_remaining=(duration/i)*(N-i)
end

if (exist(filename,'file')>0)
    known=csvread(filename);
    M=size(known,2);
    result=sortrows(vertcat(known,result));
    csvwrite(filename,result);
    
else
    csvwrite(filename,result);
end


u=csvread(filename);
v=u(:,1)
w=u(:,2)
F1_result=csvread('resonance_frequency_F1_bis.csv');
F9_result=csvread('resonance_frequency_F9.csv');
F10_result=csvread('resonance_frequency_F10.csv');
Fo_result=csvread('resonance_frequency_optimal.csv');

h_F1=F1_result(:,1);
h_F9=F9_result(:,1);
h_F10=F10_result(:,1);
h_Fo=Fo_result(:,1);

res_F1=F1_result(:,2);
res_F9=F9_result(:,2);
res_F10=F10_result(:,2);
res_Fo=Fo_result(:,2);

diff_F1=(1/hstep)*smooth(diff(smooth(res_F1,101)),151);
diff_F9=(1/hstep)*smooth(diff(smooth(res_F9,101)),151);
diff_F10=(1/hstep)*smooth(diff(smooth(res_F10,101)),151);
diff_Fo=(1/hstep)*smooth(diff(smooth(res_Fo,101)),151);

coupling_F1=smooth(diff_F1./h_F1(2:end),151);
coupling_F9=smooth(diff_F9./h_F9(2:end),151);
coupling_F10=smooth(diff_F10./h_F10(2:end),151);
coupling_Fo=smooth(diff_Fo./h_Fo(2:end),151);

% figure(1)
% hold on
% plot(freq_F1,res_F1)
% plot(freq_F9,res_F9)
% plot(freq_F10,res_F10)
% hold off
% 
% figure(2)
% hold on
% plot(freq_F1(2:end),diff_F1)
% plot(freq_F9(2:end),diff_F9)
% plot(freq_F10(2:end),diff_F10)
% hold off

h_center=0.05;
h_span=0.002;

% coupling_F1=coupling_F1((h_F1<(h_center+h_span/2)))
% coupling_F9=coupling_F9((h_F9<(h_center+h_span/2)))
% coupling_F10=coupling_F10((h_F10<(h_center+h_span/2)))
% 
% h_F1=h_F1((h_F1>(h_center-h_span/2)))
% h_F9=h_F9((h_F9>(h_center-h_span/2)))
% h_F10=h_F10((h_F10>(h_center-h_span/2)))
% pause

figure(1)
hold on
ax_Fo=plot(h_Fo((h_Fo<(h_center+h_span/2))&(h_Fo>(h_center-h_span/2))),coupling_Fo((h_Fo<(h_center+h_span/2))&(h_Fo>(h_center-h_span/2))),'-m')
ax_F1=plot(h_F1((h_F1<(h_center+h_span/2))&(h_F1>(h_center-h_span/2))),coupling_F1((h_F1<(h_center+h_span/2))&(h_F1>(h_center-h_span/2))),'-g')
ax_F9=plot(h_F9((h_F9<(h_center+h_span/2))&(h_F9>(h_center-h_span/2))),coupling_F9((h_F9<(h_center+h_span/2))&(h_F9>(h_center-h_span/2))),'-b')
ax_F10=plot(h_F10((h_F10<(h_center+h_span/2))&(h_F10>(h_center-h_span/2))),coupling_F10((h_F10<(h_center+h_span/2))&(h_F10>(h_center-h_span/2))),'-r')

legend([ax_Fo,ax_F1,ax_F9,ax_F10],'parasitic capacitance=9fF','parasitic capacitance=21.4fF','parasitic capacitance=41fF','parasitic capacitance=48fF','location','best')
title({'coupling in function of the altitude';' of the membrane for different designs'})
xlabel('height in microns')
ylabel('Coupling in microns^-1')