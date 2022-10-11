clear,close all
tic
channel_choice = [1];
quanti = 32000;
tx_amp_list = [0.1,0.15,0.2,0.25,0.35,0.5,0.65,0.8,1];

batch_size = 5000;
batch_num = 10000;

test_batch_size = 10000;
test_batch_num = 500;

M=16
%%
for tx_amp_idx=1:length(tx_amp_list)
    data_path = "./data_set_1014/train_set/pam"+M+"/tx_amp"+tx_amp_idx;
    if(~exist(data_path,'dir'))
        mkdir(char(data_path));
    end
    tx_amp = single(tx_amp_list(tx_amp_idx));
    for b_idx = 1:batch_num
        x_data = randi([0,M-1],[batch_size,1]);
        txG = real(pammod(x_data,M));

        Tx{1} = tx_amp_list(tx_amp_idx)*txG/max(txG);
        Tx{1} = round(Tx{1}*quanti);
        
        X = TxDataSort(Tx);
        trans2bin(X);
        startUDP();
        Rx = bin2receive(channel_choice);
        Rx = Rx/(-8192);

        Rx = single(Rx);
        txG = single(txG);
        
        save(data_path+"/pam16_"+b_idx+".mat","Rx","txG","tx_amp")

    end
end

%%
for tx_amp_idx=1:length(tx_amp_list)
    data_path = "./data_set_1014/test_set/pam"+M+"/tx_amp"+tx_amp_idx;
    if(~exist(data_path,'dir'))
        mkdir(char(data_path));
    end
    tx_amp = single(tx_amp_list(tx_amp_idx));
    for test_idx =1:test_batch_num
        x_data = randi([0,M-1],[test_batch_size,1]);
        txG = real(pammod(x_data,M));

        Tx{1} = tx_amp_list(tx_amp_idx)*txG/max(txG);
        Tx{1} = round(Tx{1}*quanti);
        
        X = TxDataSort(Tx);
        trans2bin(X);
        startUDP();
        Rx = bin2receive(channel_choice);
        Rx = Rx/(-8192);

        Rx = single(Rx);
        txG = single(txG);
        
        save(data_path+"/pam16_"+test_idx+".mat","Rx","txG","tx_amp")
    end
    
end

toc
