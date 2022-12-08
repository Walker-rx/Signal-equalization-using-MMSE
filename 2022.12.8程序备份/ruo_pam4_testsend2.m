% clearvars -except channel_choice bias dir_up bias_name  current dp821A  ...
%          ps upf_transmit dof_transmit filter_transmit filter_ord data_length pilot pilot_bpsk

send_length2 = length(signal_ori2);
signal_upsample = ruo_sam_rate_con(signal_ori2,filter_transmit,upf_transmit,dof_transmit);

signal_send = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp-1);
% signal_send = signal_send./norm(signal_send,2)*sqrt(length(signal_send));
signal_send = [rand([1,512]) zeros(1,10) signal_send];
Tx{1} = signal_send;

X = ruo_TxDataSort(Tx);
ruo_trans2bin(X);
ruo_startUDP();
Rx = ruo_bin2receive(channel_choice);
Rx = Rx/(-8192);

Rx2 = ruo_bin2receive(channel_choice2);
Rx2 = Rx2/(-8192);

signal_pass_channel = Rx.';
signal_pass_channel2 = Rx2.';
% signal_pass_channel = signal_send/8192;

% Rx = single(Rx);
% data_mpam = single(data_mpam);
% 
% save(data_path+"/pam"+M+"_"+"direct.mat","Rx","data_mpam")

