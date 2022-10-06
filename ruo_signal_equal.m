function x_demod = ruo_signal_equal(signal_origin,signal_current,times,coar_syn_point,pilot_length,num_of_windows,type)
    signal_ori = signal_origin(1:pilot_length); 
    signal_ori2 = upsample(signal_ori,times);
    if coar_syn_point*times-(num_of_windows*times-1) > 0
        signal_rec2 = signal_current(coar_syn_point*times-(num_of_windows*times-1):coar_syn_point*times+(num_of_windows*times)+length(signal_ori2)-1);
        [r2,~] = xcorr(signal_rec2,signal_ori2,times*2*num_of_windows-1);
        r3 = abs(r2(2*num_of_windows*times-1+1:end));
        [~,index2] = max(r3);
        index3 = index2-1+coar_syn_point*times-(num_of_windows*times-1);             
    else
        signal_rec2 = signal_current(1:coar_syn_point*times+(num_of_windows*times-1)+length(signal_ori2)-1);
        [r2,~] = xcorr(signal_rec2,signal_ori2,coar_syn_point*times+(num_of_windows*times-1)-1);
        r3 = abs(r2(coar_syn_point*times+(num_of_windows*times-1):end));
        [~,index2] = max(r3);
        index3 = index2-1+1;
    end
    fin_syn_point = index3;
    signal_received = signal_current(fin_syn_point:times:end);
%     signal_received = signal_received(1:3000);
% % % % % % % % % % % % % % %  Get H % % % % % % % % % % % % % % %
    r_phase = mod(index2,times);        % Find the phase of maximum correlation point
        if r_phase == 0
            r_phase = times;
        end
        r_for_judge = r3(r_phase:times:end);
        [~,maxlocation] = max(r_for_judge);
        if maxlocation >= 12
            channel_h_coefficient = r_for_judge(maxlocation-10:maxlocation+19);
        else
            channel_h_coefficient = r_for_judge(5:14);                          % Find 10 points with maximum correlation as h1 ... h10
        end
%         [corr_degree,~] = sort(r_for_judge,'descend');
%         channel_h_coefficient = r_for_judge(1:30);
%         channel_h_coefficient = corr_degree(1:30);  

    x_demod = 0;
    if strcmpi(type{1},'MRC')
        if maxlocation >= 12
            location_in_window = (maxlocation-10:maxlocation+19);
        else
            location_in_window = (5:14);                          
        end
    
        if coar_syn_point*times-(num_of_windows*times-1) > 0
            location_in_signal = r_phase+(location_in_window-1)*times+coar_syn_point*times-(num_of_windows*times-1)-1;
        else
            location_in_signal = r_phase+(location_in_window-1)*times;      % Find the sampling point of h1 ... h30 in received signal
        end
    
        y_phase = mod(location_in_signal,times);  % y_phase = r_phase
        if y_phase == 0
            y_phase = times;
        end
        y = signal_current(y_phase(1):times:end); % Received signal represented by symbols
        location_in_y = fix(location_in_signal/times)+1; % Find the symbol point of h1 ... h30 in received signal
    
        x_hat = 0;
        for mrc_i = 1:length(channel_h_coefficient)
            y_delay = y(location_in_y(mrc_i):location_in_y(mrc_i)+length(y)-max(location_in_y));
            x_hat = x_hat+channel_h_coefficient(mrc_i)*[y_delay,zeros(1,length(signal_origin)-length(y_delay))]; % MRC
        end
        x_hat = x_hat./norm(x_hat,2)*sqrt(length(x_hat))*sqrt(bandpower(signal_origin)); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
    
        pilot = x_hat(1:pilot_length);
        pilot = pilot./norm(pilot,2)*sqrt(length(pilot));
        pilot(pilot <= 0) = 0;
        pilot(pilot > 0) = 1;
        zero = zeros(1,10);
        signal_data = x_hat(pilot_length+11:pilot_length+11+9999);
        signal_data = signal_data./norm(signal_data,2)*sqrt(length(signal_data))*sqrt(5); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
        signal_data = (signal_data+3)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data) = 3;
        x_demod = [ pilot, zero, signal_data ];
    
    elseif strcmpi(type{1},'ZF')
% % % % % % % % % % % ZF begin % % % % % % % % % % %
%         channel_h = [fliplr(channel_h_coefficient) zeros(1,length(signal_received)-length(channel_h_coefficient)*2+1)];
%         cloumn = [channel_h(1) zeros(1,length(channel_h)-1)];
%         toep_h1 = toeplitz(cloumn,channel_h);
%         toep_h2 = [toep_h1(1:length(channel_h_coefficient)-1,...
%                         length(channel_h_coefficient):length(channel_h_coefficient)+length(channel_h_coefficient)-2),...
%                  zeros(length(channel_h_coefficient)-1,...
%                       size(toep_h1,2)-(length(channel_h_coefficient)-1))];
%         toep_h = [toep_h2;toep_h1];                                            %  H is m*n , m > n
%         tstart = tic;
%         w = toep_h'*toep_h\toep_h';                                           % equal to w = inv(toep_h'*toep_h)*toep_h';
%         toc(tstart)
%         x_hat = w*signal_received.';

        channel_h = [fliplr(channel_h_coefficient) zeros(1,length(signal_ori)-length(channel_h_coefficient))];
        cloumn = [channel_h(1) zeros(1,length(signal_received)-1)];
        toep_h = toeplitz(cloumn,channel_h);
%         channel_h = [channel_h_coefficient zeros(1,length(signal_received)-length(channel_h_coefficient))];
%         cloumn = [channel_h(1) zeros(1,length(signal_received)-1)];
%         toep_h = toeplitz(cloumn,channel_h);
        toep_h_inv = pinv(toep_h);                                           %  H is m*n , m = n
        zf_filter = fliplr(toep_h_inv(1,:));

        x_hat = conv(signal_received,zf_filter);
        x_hat = x_hat(min(length(signal_received),length(zf_filter)):end);         
%         [x_hat,~] = lteEqualizeZF(y,toep_h);
% % % % % % % % % % % ZF end % % % % % % % % % % %   
        pilot = x_hat(1:pilot_length);
        pilot = pilot./norm(pilot,2)*sqrt(length(pilot));
        pilot(pilot <= 0) = 0;
        pilot(pilot > 0) = 1;
        zero = zeros(1,10);
        signal_data = x_hat(pilot_length+11:pilot_length+11+9999);
        signal_data = signal_data./norm(signal_data,2)*sqrt(length(signal_data))*sqrt(5); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
        signal_data = (signal_data+3)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data) = 3;
        x_demod = [ pilot, zero, signal_data ];
        
    elseif strcmpi(type{1},'MMSE')      
% % % % % % % % % % % MMSE begin % % % % % % % % % % %
%         channel_h = [fliplr(channel_h_coefficient) zeros(1,length(signal_received)-length(channel_h_coefficient)*2+1)];
%         cloumn = [channel_h(1) zeros(1,length(channel_h)-1)];
%         toep_h1 = toeplitz(cloumn,channel_h);
%         toep_h2 = [toep_h1(1:length(channel_h_coefficient)-1,...
%                         length(channel_h_coefficient):length(channel_h_coefficient)+length(channel_h_coefficient)-2),...
%                  zeros(length(channel_h_coefficient)-1,...
%                       size(toep_h1,2)-(length(channel_h_coefficient)-1))];
%         toep_h = [toep_h2;toep_h1];                                            %  H is m*n , m > n
%         tstart = tic;       
        channel_h = [fliplr(channel_h_coefficient) zeros(1,length(signal_received)-length(channel_h_coefficient))];
        cloumn = [channel_h(1) zeros(1,length(signal_received)-1)];
        toep_h = toeplitz(cloumn,channel_h);
        half_mmse = (toep_h*toep_h'+type{3}/type{2}*eye(size(toep_h*toep_h',1)));
%         half_mmse = gather(half_mmse);
        half_mmse_inv = inv(half_mmse);
%         inv_half_mmse = gpuArray(inv_half_mmse);
%         w = toep_h'/(toep_h*toep_h'+type{3}/type{2}*eye(size(toep_h*toep_h',1)));
        w = toep_h'/half_mmse_inv;
%         toc(tstart)
        x_hat = w*signal_received.';
        x_hat = x_hat.';
% % % % % % % % % % % MMSE end % % % % % % % % % % %    
        pilot = x_hat(4:pilot_length+3);
        pilot = pilot./norm(pilot,2)*sqrt(length(pilot));
        pilot(pilot <= 0) = 0;
        pilot(pilot > 0) = 1;
        zero = zeros(1,10);
        signal_data = x_hat(pilot_length+14:end);
        signal_data = signal_data./norm(signal_data,2)*sqrt(length(signal_data))*sqrt(5); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
        signal_data = (signal_data+3)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data) = 3;
        x_demod = [ pilot, zero, signal_data ]; 

    elseif strcmpi(type{1},'LS')
        equal_order = type{2};
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
        pilot = pilot./norm(pilot,2)*sqrt(pilot_length)*sqrt(type{2});
        pilot(pilot <= 0) = 0;
        pilot(pilot > 0) = 1;
        zero = zeros(1,10);
        signal_data = x_hat(pilot_length+11:end);
        signal_data = signal_data./norm(signal_data,2)*sqrt(length(signal_data))*sqrt(type{3});
        signal_data = (signal_data+3)/2;
        signal_data(signal_data <= 0.5) = 0;
        signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
        signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
        signal_data(2.5 < signal_data) = 3;
        x_demod = [ pilot, zero, signal_data ];

    end
%     x_demod = x_hat;
end