% clearvars -except channel_choice bias dir_up bias_name  current dp821A  ...
%          ps upf_transmit dof_transmit filter_transmit filter_ord data_length pilot pilot_bpsk

data_path = dir_up+"train_set/pam"+M+"/bias"+bias_name+"mA/ruo"+"ruo_pam4";
if(~exist(data_path,'dir'))
    mkdir(char(data_path));
end

% data = randi([0,M-1],[1,data_length]);
% data_mpam = real(pammod(data,M));
% data_mpam_ps = bandpower(data_mpam);
% 
% unmod = [pilot,zeros(1,zero_length),data];
% signal_ori = [pilot_bpsk,zeros(1,zero_length),data_mpam];
% signal_upsample = ruo_sam_rate_con(signal_ori,filter_transmit,upf_transmit,dof_transmit);
% 
% signal_send = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp-1);
% % signal_send = signal_send./norm(signal_send,2)*sqrt(length(signal_send));
% signal_send = [rand([1,512]) zeros(1,10) signal_send];
Tx{1} = signal_send_inf;
Tx{2} = signal_send;

X_realsend = ruo_TxDataSort(Tx);
ruo_trans2bin(X_realsend);
ruo_startUDP();
Rx_realsend = ruo_bin2receive(channel_choice);
Rx_realsend = Rx_realsend/(-8192);

Rx_infsend = ruo_bin2receive(channel_choice_inf);
Rx_infsend = Rx_infsend/(-8192);

signal_pass_channel = Rx_realsend.';
signal_pass_channel_inf = Rx_infsend.';

% Rx = single(Rx);
% data_mpam = single(data_mpam);
% 
% save(data_path+"/pam"+M+"_"+"direct.mat","Rx","data_mpam")

