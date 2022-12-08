clear;
close all;

channel_choice = 1; % % 1 is direct channel
channel_choice_inf = 4;
dir_up = "./data_set_final/";
test = [1 0 1 0 1 1 1];
% bias = input('bias(mA): ')

% dp821A = visa('ni','USB0::0x1AB1::0x0E11::DP8E163250125::INSTR');
% % dp821A = visa('ni','ASRL3:INSTR');
% fopen(dp821A);
% % fprintf(dp821A,'*IDN?'); 
% % aa = fscanf(dp821A);
% % display(aa);
% fprintf(dp821A,':OUTP:OCP:VAL CH2,1.500');
% fprintf(dp821A,':OUTP:OCP:VAL? CH2');
% aa = fscanf(dp821A);
% display(aa);
% fprintf(dp821A,':OUTP:OCP CH2,ON');
% fprintf(dp821A,':OUTP CH2,ON');
% fprintf(dp821A,':APPL CH1,24,1.000');
% fprintf(dp821A,':OUTP CH1,ON');

bias_name = 780;
pause(0.5);

M = 4;
data_length = 10000;
zero_length = 3000;
ls_order = 50;

origin_rate = 25e6; 
bw = origin_rate/2;         % baseband bandwidth
f_rate = 160e6;
d_rate = 150e6;

times = d_rate/origin_rate;
num_of_windows = 100;

% pilot = pilot_gen([0 0 0 1 1 1 1]);
% pilot = pilot_gen([1 1 1 0 0 0 1 1 1]);
pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 0 1]);    %  2^11 = 2048
% pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 1 0 1 1]);   %  2^13 = 8192
pilot_ps = bandpower(pilot);
pilot_length = length(pilot);
pilot_bpsk = pilot*2-1;

filter_ord = 1000;     % filter order used in function sam_rate_con
rp = 0.00057565;      
rst = 1e-4;       % filter parameter used in function sam_rate_con

llcm_transmit = lcm(origin_rate,f_rate);
upf_transmit = llcm_transmit/origin_rate;  % upsampling parameters
dof_transmit = llcm_transmit/f_rate;
ups_rate_transmit = origin_rate*upf_transmit;
filter_transmit = firceqrip(filter_ord,bw/(ups_rate_transmit/2),[rp rst],'passedge'); % Impulse function
filter_transmit = filter_transmit./norm(filter_transmit,2)*sqrt(bw/ups_rate_transmit*2); % filter used in function sam_rate_con in transmit part

llcm_receive = lcm(f_rate,d_rate);
upf_receive = llcm_receive/f_rate;  % upsampling parameters
dof_receive = llcm_receive/d_rate;
ups_rate_receive = f_rate*upf_receive;
filter_receive = firceqrip(filter_ord,bw/(ups_rate_receive/2),[rp rst],'passedge'); % Impulse function
filter_receive = filter_receive./norm(filter_receive,2)*sqrt(bw/ups_rate_receive*2); % filter used in function sam_rate_con in receive part

% origin_rate = gpuArray(double(origin_rate));
% f_rate = gpuArray(double(f_rate));
% d_rate = gpuArray(double(d_rate));
% num_of_windows = gpuArray(double(num_of_windows));
% times = gpuArray(double(times));
% pilot_length = gpuArray(double(pilot_length));
% delay = gpuArray(double(delay));
% filter_transmit = gpuArray(double(filter_transmit));
% filter_receive = gpuArray(double(filter_receive));
% upf_transmit = gpuArray(double(upf_transmit));
% dof_transmit = gpuArray(double(dof_transmit));
% upf_receive = gpuArray(double(upf_receive));
% dof_receive = gpuArray(double(dof_receive));
% h_channel = gpuArray(double(h_channel));
% h_channel_delay = gpuArray(double(h_channel_delay));
   
amp = 50;
looptime = 21;

ruo_load_data

data = data_demod_nofalse;
data_mpam = real(pammod(data,M));
data_mpam_ps = bandpower(data_mpam);

unmod = [pilot,zeros(1,zero_length),data];

%% Locating the transmission error point

signal_demod_ls_inf = ruo_signal_equal(signal_ori,signal_received_inf,times,coar_syn_point_inf,pilot_length,zero_length,num_of_windows,{'ls',ls_order});

data_demod_ls_inf = signal_demod_ls_inf(pilot_length+zero_length+1:end);

if length(signal_demod_ls_inf) < length(unmod)
    compare_length_inf = length(signal_demod_ls_inf);
    error_location_inf = find(signal_demod_ls_inf ~= unmod(1:compare_length_inf));
else
    compare_length_inf = length(unmod);
    error_location_inf = find(signal_demod_ls_inf(1:compare_length_inf) ~= unmod);
end
txerror_num_loop = length(error_location_inf) ;
%% Calculating the SER without correction
signal_received_channel1_send1 = signal_received_notcorrect;
[length_loop_channel1_send1 , ps_loop_channel1_send1 , pn_loop_channel1_send1 , ...
    errornum_ls_loop_channel1_send1 , error_location_loop_channel1_send1 , data_demod_ls_channel1_send1] ...
    = ruo_calculate_ser( data , signal_ori , signal_received_channel1_send1 , pilot_length , zero_length , data_length ...
    , coar_syn_point , fin_syn_point , times , num_of_windows , ls_order);

signal_received_channel1_send2 = signal_received_forcorrect;
[length_loop_channel1_send2 , ps_loop_channel1_send2 , pn_loop_channel1_send2 , ...
    errornum_ls_loop_channel1_send2 , error_location_loop_channel1_send2 , data_demod_ls_channel1_send2] ...
    = ruo_calculate_ser( data , signal_ori , signal_received_channel1_send2 , pilot_length , zero_length , data_length ...
    , coar_syn_point_forcorrect , fin_syn_point_forcorrect , times , num_of_windows , ls_order);

%% Replacing the transmission error point
signal_received_tmp = signal_received_channel1_send1;
signal_received_aftercorrect = signal_received_channel1_send1;
fin_syn_point_forcorrect = ruo_signal_syn_recorrect(signal_received_aftercorrect(fin_syn_point:end),signal_received_channel1_send2,fin_syn_point_forcorrect);
replace_loc = error_location_inf*times + fin_syn_point - 1;
replace_loc_correct = error_location_inf*times + fin_syn_point_forcorrect - 1;
replace_length = 8*times;
for i = 1:txerror_num_loop
    signal_received_aftercorrect(replace_loc(i)-replace_length : replace_loc(i)+replace_length) = signal_received_channel1_send2(replace_loc_correct(i)-replace_length : replace_loc_correct(i)+replace_length);
end
%% Calculating the SER after correction
[length_loop_aftercorrect , ps_loop_aftercorrect , pn_loop_aftercorrect , ...
    errornum_ls_loop_aftercorrect , error_location_aftercorrect , data_demod_ls_aftercorrect] ...
    = ruo_calculate_ser( data , signal_ori , signal_received_aftercorrect , pilot_length , zero_length , data_length ...
    , coar_syn_point, fin_syn_point , times , num_of_windows , ls_order);
%%
fprintf('%d %d %d \n',errornum_ls_loop_channel1_send1,txerror_num_loop,errornum_ls_loop_aftercorrect);


