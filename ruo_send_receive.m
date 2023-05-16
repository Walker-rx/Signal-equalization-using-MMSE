clear;
close all;

%% Initial parameters
channel_choice = 1; 
channel_choice_inf = 3;

M = 2;
data_length_initial = 10000;
zero_length = 3000;
zero_length_forsyn = 200;
ls_order = 50;
num_of_windows = 100;

times = 6;
origin_rate_tmp = 10e6;
f_rate = 160e6;
d_rate_tmp = origin_rate_tmp*times;

filter_order = 2000;     % filter order used in function sam_rate_con
rp = 0.00057565;      
rst = 1e-4;       % filter parameter used in function sam_rate_con

%% Generate filter
[origin_rate , d_rate , upf_transmit , dof_transmit , filter_transmit , upf_receive , dof_receive , filter_receive] = ruo_filter_gen(origin_rate_tmp , f_rate , times , filter_order , rp , rst);
rate_change = 100*abs(origin_rate-origin_rate_tmp)/origin_rate_tmp;
fprintf("rate's change = %.8f%% \n",rate_change);
ups_time = upf_transmit/dof_transmit;
% data_length = round(data_length_initial*8/ups_time);
data_length = data_length_initial;

filter_10m = firceqrip(filter_order,0.5e6/(origin_rate/2),[rp rst],'high');
%% Generate pilot
pilot_ini = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 0 1]);    %  2^11 = 2048
pilot_bpsk_ini = pilot_ini*2-1;
pilot_bpsk_tmp = conv(pilot_bpsk_ini,filter_10m);
pilot_bpsk_tmp = pilot_bpsk_tmp((length(filter_10m)+1)/2 : length(pilot_bpsk_tmp)-(length(filter_10m)-1)/2);
pilot_bpsk_tmp(pilot_bpsk_tmp <= 0) = -1;
pilot_bpsk_tmp(pilot_bpsk_tmp > 0) = 1;
pilot_bpsk = pilot_bpsk_tmp;
pilot = ruo_pamdemod(pilot_bpsk,2);
pilot_length = length(pilot);

%% Signal amplitude range
t = datetime('now');

% amp_begin = 0.0015 , amp_end = 1 , amp_step = 0.03994
% amp_begin = 0.1613;
amp_step = 0.15976;
amp_begin = 1;
amp_end = 1;
amp_step = 0.15976;
% amp_begin = 0.48082;
% amp_end = 0.48082;
% amp_step = 1;
amp_inf = 0.6;
total_num = (amp_end-amp_begin)/amp_step+1;

fprintf('bpsk data \n');
%% Current range
% dp821A = visa('ni','USB0::0x1AB1::0x0E11::DP8D163650047::INSTR');
dp821A = visa('ni','TCPIP::121.248.13.210::INSTR');


bias_begin = 50;
bias_step = 40;
bias_end = 850;
save_num = 400;

%% 
save_amp_num = 0;
for amp = amp_begin:amp_step:amp_end
    fopen(dp821A);
    save_amp_num = save_amp_num + 1;
    save_path_ini = "dnn_data/"+t.Year+"."+t.Month+"."+t.Day+"/"+origin_rate_tmp/1e6+"M/amp"+amp;
    for current = bias_begin : bias_step : bias_end
        fprintf(dp821A,strcat(':APPL CH1,10,',num2str(current/1000)));
        fprintf(dp821A,':OUTP CH1,ON');
        pause(3);
        fprintf(dp821A,':MEAS:CURR? CH1');
        bias = str2num(fscanf(dp821A));
        pause(3);
        
        save_path = save_path_ini+"/bias"+current/1000;
        if(~exist(save_path,'dir'))
            mkdir(char(save_path));
        end
        %%
        ruo_gen_light_data
        
        %%
        if current == bias_begin
            save_bias = fopen(save_path_ini+"/save_bias.txt",'w');
            
            save_parameter = fopen(save_path_ini+"/save_parameter.txt",'w');
            fprintf(save_parameter,' amp begin = %f , amp end = %f , amp step = %f \n',amp_begin,amp_end,amp_step);
            fprintf(save_parameter,' bias begin = %f , bias end = %f , bias step = %f \n',bias_begin,bias_end,bias_step);
            fclose(save_parameter);
        else
            save_bias = fopen(save_path_ini+"/save_bias.txt",'a');
        end
        fprintf(save_bias,'%f \n',bias);
        fclose(save_bias);
        
        pause(10);
        fprintf(dp821A,':OUTP CH1,OFF');
        pause(3);
        
    end
    fclose(dp821A);
end
%%





