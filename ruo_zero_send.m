% clearvars -except channel_choice bias dir_up bias_name current dp821A data_length 

zero = zeros(1,data_length);
Tx{1} = zero;

X = ruo_TxDataSort(Tx);
ruo_trans2bin(X);
ruo_startUDP();
Rx = ruo_bin2receive(channel_choice);
Rx = Rx/(-8192);
% Rx = Rx/(-1);

noise = Rx.';




