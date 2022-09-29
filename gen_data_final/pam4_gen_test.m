clearvars -except channel_choice quanti bias dir_up bias_name  current dp821A

seq_length = 3000;
M = 4;

data_path = dir_up+"train_set/pam"+M+"/bias"+bias_name+"mA/test"+"send_receive_test";
if(~exist(data_path,'dir'))
    mkdir(char(data_path));
end

x_data = randi([0,M-1],[seq_length,1]);
txG = real(pammod(x_data,M));

Tx_tem = 0.5*txG/max(txG);
Tx_tem = round(Tx_tem*quanti);
Tx{1} = conv(Tx_tem,filter_transmit);

X = TxDataSort(Tx);
trans2bin_test(X);
startUDP_test();
Rx = bin2receive_test(channel_choice);
Rx = Rx/(-8192);

Rx = single(Rx);
txG = single(txG);

save(data_path+"/pam"+M+"_"+"test.mat","Rx","txG","0.5","bias")

