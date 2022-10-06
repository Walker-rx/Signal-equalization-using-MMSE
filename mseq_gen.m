clearvars -except channel_choice quanti bias dir_up bias_name current dp821A
tic
load("msequence.mat")

tx_amp_list = [0.1,0.15,0.25,0.35,0.5,0.65,0.8,1];

for tx_amp_idx=1:length(tx_amp_list)
    data_path = dir_up+"train_set/mseq/bias"+bias_name+"mA/tx_amp"+tx_amp_idx;
    if(~exist(data_path,'dir'))
        mkdir(char(data_path));
    end
    tx_amp = single(tx_amp_list(tx_amp_idx));
    txG = mseq18;
    Tx{1} = tx_amp_list(tx_amp_idx)*txG/max(txG);
    Tx{1} = round(Tx{1}*quanti);

    X = TxDataSort(Tx);
    trans2bin(X);
    startUDP();
    Rx = bin2receive(channel_choice);
    Rx = Rx/(-8192);

    Rx = single(Rx);
    txG = single(txG);
    
    save(data_path+"/mseq18.mat","Rx","txG","tx_amp","bias")
end

for tx_amp_idx=1:length(tx_amp_list)
    data_path = dir_up+"train_set/mseq/bias"+bias_name+"mA/tx_amp"+tx_amp_idx;
    if(~exist(data_path,'dir'))
        mkdir(char(data_path));
    end
    tx_amp = single(tx_amp_list(tx_amp_idx));
    txG = mseq16;
    Tx{1} = tx_amp_list(tx_amp_idx)*txG/max(txG);
    Tx{1} = round(Tx{1}*quanti);

    X = TxDataSort(Tx);
    trans2bin(X);
    startUDP();
    Rx = bin2receive(channel_choice);
    Rx = Rx/(-8192);

    Rx = single(Rx);
    txG = single(txG);
    
    save(data_path+"/mseq16.mat","Rx","txG","tx_amp","bias")
end

for tx_amp_idx=1:length(tx_amp_list)
    data_path = dir_up+"train_set/mseq/bias"+bias_name+"mA/tx_amp"+tx_amp_idx;
    if(~exist(data_path,'dir'))
        mkdir(char(data_path));
    end
    tx_amp = single(tx_amp_list(tx_amp_idx));
    txG = mseq17;
    Tx{1} = tx_amp_list(tx_amp_idx)*txG/max(txG);
    Tx{1} = round(Tx{1}*quanti);

    X = TxDataSort(Tx);
    trans2bin(X);
    startUDP();
    Rx = bin2receive(channel_choice);
    Rx = Rx/(-8192);

    Rx = single(Rx);
    txG = single(txG);
    
    save(data_path+"/mseq17.mat","Rx","txG","tx_amp","bias")
end

for tx_amp_idx=1:length(tx_amp_list)
    data_path = dir_up+"train_set/mseq/bias"+bias_name+"mA/tx_amp"+tx_amp_idx;
    if(~exist(data_path,'dir'))
        mkdir(char(data_path));
    end
    tx_amp = single(tx_amp_list(tx_amp_idx));
    txG = mseq14;
    Tx{1} = tx_amp_list(tx_amp_idx)*txG/max(txG);
    Tx{1} = round(Tx{1}*quanti);

    X = TxDataSort(Tx);
    trans2bin(X);
    startUDP();
    Rx = bin2receive(channel_choice);
    Rx = Rx/(-8192);

    Rx = single(Rx);
    txG = single(txG);
    
    save(data_path+"/mseq14.mat","Rx","txG","tx_amp","bias")
end

disp("mseq completed.")
toc
