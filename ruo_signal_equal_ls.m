function x_demod = ruo_signal_equal_ls(signal_origin,signal_current,times,fin_syn_point,pilot_length,zero_length,equal_order,M)
    w = zeros(equal_order,times);
    U = zeros(equal_order,pilot_length-equal_order+1);
    Err = zeros(1,times);
    %         start_point = fin_syn_point - 9*times;
    headwindow = equal_order-(fix(equal_order/2)+1);

    if headwindow*times > fin_syn_point-1
        signal_addzero = [ zeros(1,headwindow*times-(fin_syn_point-1)) , signal_current ];
    else
        start_point = fin_syn_point - headwindow*times;
        signal_addzero = signal_current(start_point:end);
    end

    for i = 0:(times-1)
        pilot_received = downsample(signal_addzero(1:pilot_length*times),times,i);
        %             pilot_received = downsample(signal_current(start_point-1 + (1:pilot_length*times)),times,i);
        for j = equal_order:pilot_length
            U(:,j-equal_order+1) = pilot_received(j:-1:j-equal_order+1)';
        end
        pilot_ori = signal_origin(1:size(U,2))';
        w(:,i+1) = (U*U')\U*pilot_ori;
        error = U'*w(:,i+1)-pilot_ori;
        Err(i+1) = norm(error,2);
    end
    [~,phase] = min(Err);
    %         signal_received_ls = downsample(signal_current(start_point:end),times,phase-1);
    signal_received_ls = downsample(signal_addzero,times,phase-1);
    W = w(:,phase);
    x_hat = conv(signal_received_ls,W,'valid');

    pilot = x_hat(1:pilot_length);
    pilot(pilot <= 0) = 0;
    pilot(pilot > 0) = 1;
%     pilot = ruo_pamdemod(pilot,M);
    zero = zeros(1,zero_length);
    signal_data_tmp = x_hat(pilot_length+zero_length+1:end);
    %         signal_data(signal_data <= 0) = 0;
    %         signal_data(signal_data > 0) = 1;
    signal_data = ruo_pamdemod(signal_data_tmp,M);
    x_demod = [ pilot, zero, signal_data ];
end
