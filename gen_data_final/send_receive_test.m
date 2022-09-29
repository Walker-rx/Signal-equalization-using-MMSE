clear,close all

channel_choice = [4];
quanti = 32000;
dir_up = "./data_set_final/";
% bias = input('bias(mA): ')

dp821A = visa('ni','USB0::0x1AB1::0x0E11::DP8E163250125::INSTR');
% dp821A = visa('ni','ASRL3:INSTR');
fopen(dp821A);
% fprintf(dp821A,'*IDN?');
% aa = fscanf(dp821A);
% display(aa);
fprintf(dp821A,':OUTP:OCP:VAL CH2,1.500');
fprintf(dp821A,':OUTP:OCP:VAL? CH2');
aa = fscanf(dp821A);
display(aa);
fprintf(dp821A,':OUTP:OCP CH2,ON');
fprintf(dp821A,':OUTP CH2,ON');
fprintf(dp821A,':APPL CH1,24,1.000');
fprintf(dp821A,':OUTP CH1,ON');

origin_rate = 25e6; 
bw = origin_rate/2;         % baseband bandwidth
f_rate = 160e6;
d_rate = 150e6;

filter_ord = 1000;     % filter order used in function sam_rate_con
rp = 0.00057565;      
rst = 1e-4;       % filter parameter used in function sam_rate_con

llcm_transmit = lcm(origin_rate,f_rate);
upf_transmit = llcm_transmit/origin_rate;  % upsampling parameters
ups_rate_transmit = origin_rate*upf_transmit;
filter_transmit = firceqrip(filter_ord,bw/(ups_rate_transmit/2),[rp rst],'passedge'); % Impulse function
filter_transmit = filter_transmit./norm(filter_transmit,2)*sqrt(bw/ups_rate_transmit*2); % filter used in function sam_rate_con in transmit part

for current=780:40:820
    bias_name = current
    fprintf(dp821A,strcat(':APPL CH2,8,',num2str(current/1000)));
    pause(1.5);
    fprintf(dp821A,':MEAS:CURR? CH2');
    aa = fscanf(dp821A);
    bias = str2num(aa)
    pause(0.5);
    %%
    pam4_gen_test
    
end
fprintf(dp821A,':OUTP CH2,OFF');
fprintf(dp821A,':OUTP CH1,OFF');
fclose(dp821A);