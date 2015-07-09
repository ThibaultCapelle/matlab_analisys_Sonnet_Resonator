filename='resonance_A1.csv';
file = fopen(filename,'r');
fgetl(file);
freq=[];
S21=[];
i=0;
while ~feof(file)
    i=i+1;
    line=fgetl(file);
    line_s=strsplit(line);
    f=line_s(1);
    f=str2num(f{1});
    s=line_s(2);
    s=str2num(s{1});
    freq=vertcat(freq,f);
    S21=vertcat(S21,s);
end
N=i;

result=horzcat(freq,real(S21),imag(S21));
csvwrite('resonance_A1_formatted.csv',result)
fclose(file);