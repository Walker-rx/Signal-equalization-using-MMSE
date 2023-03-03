function data_demod = ruo_pamdemod(signal_origin,M) 
    if M == 2
        signal_data = signal_origin;
        signal_data(signal_data <= 0) = 0;
        signal_data(signal_data > 0) = 1;       
    elseif M == 4
        signal_data = (signal_origin+M-1)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data) = 3;
    elseif M == 5
        signal_data = (signal_origin+M-1)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data & signal_data <= 3.5) = 3;
        signal_data(3.5 < signal_data) = 4;
    elseif M == 6
        signal_data = (signal_origin+M-1)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data & signal_data <= 3.5) = 3;
        signal_data(3.5 < signal_data & signal_data <= 4.5) = 4;
        signal_data(4.5 < signal_data) = 5;
    elseif M == 7
        signal_data = (signal_origin+M-1)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data & signal_data <= 3.5) = 3;
        signal_data(3.5 < signal_data & signal_data <= 4.5) = 4;
        signal_data(4.5 < signal_data & signal_data <= 5.5) = 5;
        signal_data(5.5 < signal_data) = 6;
    elseif M == 8
        signal_data = (signal_origin+M-1)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data & signal_data <= 3.5) = 3;
        signal_data(3.5 < signal_data & signal_data <= 4.5) = 4;
        signal_data(4.5 < signal_data & signal_data <= 5.5) = 5;
        signal_data(5.5 < signal_data & signal_data <= 6.5) = 6;
        signal_data(6.5 < signal_data) = 7;
    elseif M == 16
        signal_data = (signal_origin+M-1)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data & signal_data <= 3.5) = 3;
        signal_data(3.5 < signal_data & signal_data <= 4.5) = 4;
        signal_data(4.5 < signal_data & signal_data <= 5.5) = 5;
        signal_data(5.5 < signal_data & signal_data <= 6.5) = 6;       
        signal_data(6.5 < signal_data & signal_data <= 7.5) = 7;
        signal_data(7.5 < signal_data & signal_data <= 8.5) = 8;
        signal_data(8.5 < signal_data & signal_data <= 9.5) = 9;
        signal_data(9.5 < signal_data & signal_data <= 10.5) = 10;
        signal_data(10.5 < signal_data & signal_data <= 11.5) = 11;
        signal_data(11.5 < signal_data & signal_data <= 12.5) = 12;       
        signal_data(12.5 < signal_data & signal_data <= 13.5) = 13;
        signal_data(13.5 < signal_data & signal_data <= 14.5) = 14;
        signal_data(14.5 < signal_data) = 15;
    end
    data_demod = signal_data;
end