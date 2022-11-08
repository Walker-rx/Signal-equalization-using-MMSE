% clearvars -except channel_choice bias dir_up bias_name  current dp821A  ...
%          ps upf_transmit dof_transmit filter_transmit filter_ord data_length pilot pilot_bpsk

% namestr_ori = ['signal_ori_save_' num2str(mat_location)];
% signal_ori = eval(namestr_ori);
data_mpam = signal_ori(pilot_length+zero_length+1:end);
data = (data_mpam+3)/2;
signal_upsample = ruo_sam_rate_con(signal_ori,filter_transmit,upf_transmit,dof_transmit);

signal_send = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp-1);
% signal_send = signal_send./norm(signal_send,2)*sqrt(length(signal_send));
signal_send = [rand([1,512]) zeros(1,10) signal_send];
Tx{1} = signal_send;

X = ruo_TxDataSort(Tx);
ruo_trans2bin(X);
ruo_startUDP();
Rx = ruo_bin2receive(channel_choice);
Rx = Rx/(-8192);
% Rx = Rx/(-1);
signal_pass_channel = Rx.';

% Rx = single(Rx);
% data_mpam = single(data_mpam);
% 
% save(data_path+"/pam"+M+"_"+"direct.mat","Rx","data_mpam")

