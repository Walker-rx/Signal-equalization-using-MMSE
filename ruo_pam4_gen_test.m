% clearvars -except channel_choice quanti bias dir_up bias_name  current dp821A filter_transmit
rng(0);
seq_length = 10000;
M = 4;

data_path = dir_up+"train_set/pam"+M+"/bias"+bias_name+"mA/test"+"send_receive_test";
if(~exist(data_path,'dir'))
    mkdir(char(data_path));
end

x_data = randi([0,M-1],[1,seq_length]);
txG = real(pammod(x_data,M));

Tx_tem = 0.5*txG/max(txG);
Tx_tem = round(Tx_tem*quanti);
Tx_tem = upsample(Tx_tem,6);
Tx{1} = conv(Tx_tem,filter_transmit);
Tx{1} = Tx{1}((length(filter_transmit)+1)/2:length(Tx{1})-(length(filter_transmit)-1)/2);

X = ruo_TxDataSort(Tx);
ruo_trans2bin(X);
ruo_startUDP();
Rx = ruo_bin2receive(channel_choice);
Rx = Rx/(-8192);

Rx = single(Rx);
txG = single(txG);

% save(data_path+"/pam"+M+"_"+"test.mat","Rx","txG","0.5","bias")

