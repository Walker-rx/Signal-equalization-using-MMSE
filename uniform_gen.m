clearvars -except channel_choice quanti bias dir_up bias_name  current dp821A
tic


batch_size = 3000;
batch_num = 2000;

test_batch_size = 3000;
test_batch_num = 500;

tx_amp_list = [0.1,0.15,0.25,0.35,0.5,0.65,0.8,1];
%%
data_path = dir_up+"train_set/uniform/bias"+bias_name+"mA/";
if(~exist(data_path,'dir'))
    mkdir(char(data_path));
end
for b_idx = 1:batch_num
    txG = 2*rand([batch_size,1])-1;
    
    Tx{1} = txG;
    Tx{1} = round(Tx{1}*quanti);
    
    X = TxDataSort(Tx);
    trans2bin(X);
    startUDP();
    Rx = bin2receive(channel_choice);
    Rx = Rx/(-8192);
    
    Rx = single(Rx);
    txG = single(txG);
    save(data_path+"uni_"+b_idx+".mat","Rx","txG","bias")
    
end

%%

for tx_amp_idx=1:length(tx_amp_list)
    data_path = dir_up+"test_set/uniform/bias"+bias_name+"mA/tx_amp"+tx_amp_idx;
    if(~exist(data_path,'dir'))
        mkdir(char(data_path));
    end
    for test_idx =1:test_batch_num
        txG = 2*rand([test_batch_size,1])-1;
        
        Tx{1} = txG*tx_amp_list(tx_amp_idx);
        Tx{1} = round(Tx{1}*quanti);

        X = TxDataSort(Tx);
        trans2bin(X);
        startUDP();
        Rx = bin2receive(channel_choice);
        Rx = Rx/(-8192);

        Rx = single(Rx);
        txG = single(txG);
        tx_amp = single(tx_amp_list(tx_amp_idx));
        save(data_path+"/uni_"+test_idx+".mat","Rx","txG","tx_amp","bias")
    end
end
disp("uniform completed.")
toc