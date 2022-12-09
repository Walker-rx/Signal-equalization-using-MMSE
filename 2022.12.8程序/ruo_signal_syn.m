%%%%%%%%%%%%%%%%%%%%% Signal synchronization %%%%%%%%%%%%%%%%%%%%%%%%%
 
function [fin_syn_point,coar_syn_point] = ruo_signal_syn(origin_sam_rate,cur_sam_rate,origin_signal,current_signal,num_of_windows) 
%signal's bandwidth; current sampling rate; origin signal;signal after
%sampling; 8., 
    if(~exist('num_of_windows','var'))
        num_of_windows = 10;                % number of symbols moved left and right after rough sampling
    end

    times = cur_sam_rate/origin_sam_rate;
    maxlag = 500*times;
    % signal_ori = fliplr(sys_mpam(1:256));
%     signal_ori = origin_signal(1:times*50); 
    signal_ori = origin_signal(1:times*100);
    
%     curr_judge = current_signal(1:1000*times);
    curr_judge = current_signal(1:4000*times);
    energy = zeros(1,times);
    for downsample_phase = 0:times-1
        sampled_curr_judge = downsample(curr_judge,times,downsample_phase);
        energy(downsample_phase+1) = norm(sampled_curr_judge,2);
    end
    [~,sampling_phase] = max(energy);

    signal_resam = current_signal(sampling_phase:times:end);            % take a sampling point from each symbol of the signal
    if times*3000 <= length(signal_resam)
        signal_rec = signal_resam(1:times*3000);
    else
        signal_rec = signal_resam(1:end);
    end
%     signal_rec = signal_resam(1:4000);
    [r0,~] = xcorr(signal_rec,signal_ori,maxlag);    
    r1 = abs(r0(maxlag:end));
    [~,index0] = max(r1);                            % Find the point with the largest correlation coefficient
    if index0  ==  1
        coar_syn_point = 1;
    else
        coar_syn_point = index0-1;                         % Coarse synchronization point
    end

%     subplot(2,1,1);
%     plot(offset,r0);
%     xlabel('偏移量m');ylabel('相关系数');
%     title(sprintf("粗同�?%d",coar_syn_point));

    signal_ori2 = upsample(signal_ori,times);
    if coar_syn_point*times-(num_of_windows*times-1) > 0
        signal_rec2 = current_signal(coar_syn_point*times-(num_of_windows*times-1):coar_syn_point*times+(num_of_windows*times-1)+length(signal_ori2)-1);
        [r2,offset2] = xcorr(signal_rec2,signal_ori2,times*2*num_of_windows-2);
        r3 = abs(r2(2*num_of_windows*times-2+1:end));
        [~,index2] = max(r3);
        index3 = index2-1+coar_syn_point*times-(num_of_windows*times-1);
        fin_syn_point = index3;                                  % Fine synchronization

%         plot_offset = offset2-(num_of_windows*times-1);       % Move the figure to the origin point
%         fine_syn_offset = index2-(num_of_windows*times-1)-1;       
    else
        signal_rec2 = current_signal(1:coar_syn_point*times+(num_of_windows*times-1)+length(signal_ori2)-1);
        [r2,offset2] = xcorr(signal_rec2,signal_ori2,coar_syn_point*times+(num_of_windows*times-1)-1);
        r3 = abs(r2(coar_syn_point*times+(num_of_windows*times-1):end));
        [~,index2] = max(r3);
        index3 = index2-1+1;
        fin_syn_point = index3;   

%         plot_offset = offset2-(coar_syn_point*times-1);
%         fine_syn_offset = index2-(coar_syn_point*times-1)-1; 
    end
%     subplot(2,1,2);
%     plot(plot_offset,r2);
%     xlabel('偏移量m');ylabel('相关系数');
%     title([strcat("精同�?",num2str(fin_syn_point),',',"偏移�?",num2str(fine_syn_offset))]);
end


% function [syn_signal,coar_syn_point] = signal_syn(signal_bw,cur_sam_rate,origin_signal,current_signal,ps,sigma2,varargin) 
% %signal's bandwidth; current sampling rate; origin signal;signal after sampling; 
% 
% %     if(~exist('num_of_windows','var'))
% %         num_of_windows = 2;
% %     end
%     p = inputParser;
%     p.addOptional('num_of_windows',10);     % number of symbols moved left and right after rough sampling
%     p.addParameter('demodulation','False',@(x)any(validatestring(x,{'MRC','ZF','MMSE','False'})));
%     p.parse(varargin{:});
% 
%     times = cur_sam_rate/signal_bw;
%     maxlag = 50*times;
%     % signal_ori = fliplr(sys_mpam(1:256));
%     signal_ori = origin_signal(1:times*50); 
% 
%     curr_judge = current_signal(1:1000*times);
%     energy = zeros(1,times);
%     for downsample_phase = 0:times-1
%         sampled_curr_judge = downsample(curr_judge,times,downsample_phase);
%         energy(downsample_phase+1) = norm(sampled_curr_judge,2);
%     end
% 
%     [~,sampling_phase] = max(energy);
% 
%     s_resam = current_signal(sampling_phase:times:end);            % take a sampling point from each symbol of the signal
%     signal_rec = s_resam(1:times*300);
% %     disp(phase);
% 
%     [r0,offset] = xcorr(signal_rec,signal_ori,maxlag);    
%     r1 = abs(r0(maxlag+1:end));
%     [~,index0] = max(r1);                            % Find the point with the largest correlation coefficient
% %     index1 = index0-(maxlag+1)+1;
%     coar_syn_point = index0;                                 % Coarse synchronization
% 
%     subplot(2,1,1);
%     plot(offset,r0);
%     xlabel('偏移量m');ylabel('相关系数');
%     title(sprintf("粗同�?%d",coar_syn_point));
% 
%     signal_ori2 = upsample(signal_ori,times);
%     if coar_syn_point*times-(p.Results.num_of_windows*times-1) > 0
%         signal_rec2 = current_signal(coar_syn_point*times-(p.Results.num_of_windows*times-1):coar_syn_point*times+(p.Results.num_of_windows*times-1)+length(signal_ori2)-1);
%         [r2,offset2] = xcorr(signal_rec2,signal_ori2,times*2*p.Results.num_of_windows-2);
%         r3 = abs(r2(2*p.Results.num_of_windows*times-2+1:end));
%         [~,index2] = max(r3);
%         index3 = index2-1+coar_syn_point*times-(p.Results.num_of_windows*times-1);
%         syn_signal = index3;                                  % Fine synchronization
% 
%         plot_offset = offset2-(p.Results.num_of_windows*times-1);       % Move the figure to the origin point
%         fine_syn_offset = index2-(p.Results.num_of_windows*times-1)-1;       
%     else
%         signal_rec2 = current_signal(1:coar_syn_point*times+(p.Results.num_of_windows*times-1)+length(signal_ori2)-1);
%         [r2,offset2] = xcorr(signal_rec2,signal_ori2,coar_syn_point*times+(p.Results.num_of_windows*times-1)-1);
%         r3 = abs(r2(coar_syn_point*times+(p.Results.num_of_windows*times-1):end));
%         [~,index2] = max(r3);
%         index3 = index2-1+1;
%         syn_signal = index3;   
% 
%         plot_offset = offset2-(coar_syn_point*times-1);
%         fine_syn_offset = index2-(coar_syn_point*times-1)-1; 
%     end
%     subplot(2,1,2);
%     plot(plot_offset,r2);
%     xlabel('偏移量m');ylabel('相关系数');
% %     title(sprintf("精同�?%d",syn_signal),sprintf("偏移�?%d",index2-1));
%     title([strcat("精同�?",num2str(syn_signal),',',"偏移�?",num2str(fine_syn_offset))]);
% 
%     x_demod = 0;
%     if strcmpi(p.Results.demodulation,'MRC')
%         r_phase = mod(index2,times);          % Find the phase of maximum correlation point
%         r_for_judge = r3(r_phase:times:end);
%         [corr_degree,location] = sort(r_for_judge,'descend');
%         channel_h_coefficient = corr_degree(1:30);
%         location_in_window = location(1:30);  % Find 30 points with maximum correlation as h1 ... h30
% 
%         if coar_syn_point*times-(p.Results.num_of_windows*times-1) > 0
%             location_in_signal = r_phase+(location_in_window-1)*times+coar_syn_point*times-(p.Results.num_of_windows*times-1)-1;
%         else
%             location_in_signal = r_phase+(location_in_window-1)*times;      % Find the sampling point of h1 ... h30 in received signal
%         end
% 
%         y_phase = mod(location_in_signal,times);  % y_phase = r_phase
%         y = current_signal(y_phase(1):times:end); % Received signal represented by symbols
%         location_in_y = fix(location_in_signal/times)+1; % Find the symbol point of h1 ... h30 in received signal
%                 
%         x_hat = 0;
%         for mrc_i = 1:length(channel_h_coefficient)
%             x_hat = x_hat+channel_h_coefficient(mrc_i)*y(location_in_y(mrc_i):location_in_y(mrc_i)+length(y)-max(location_in_y)); % MRC
%         end
%         x_hat = x_hat./norm(x_hat,2)*sqrt(length(x_hat))*sqrt(5); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
%                  
%         x_demod = (x_hat+3)/2;
%         x_demod(x_demod <= 0.5) = 0;
%         x_demod(0.5 < x_demod & x_demod <= 1.5) = 1;
%         x_demod(1.5 < x_demod & x_demod <= 2.5) = 2;
%         x_demod(2.5 < x_demod) = 3;
% 
%     elseif strcmpi(p.Results.demodulation,'ZF')
%         r_phase = mod(index2,times);        % Find the phase of maximum correlation point
%         r_for_judge = r3(r_phase:times:end); 
%         [corr_degree,~] = sort(r_for_judge,'descend');
%         channel_h_coefficient = corr_degree(1:30);  % Find 30 points with maximum correlation as h1 ... h30
%         y = current_signal(index3:times:index3+times*(length(signal_ori)-1)); % Received signal represented by symbols.
%         % % % % % % % % % % % ZF begin % % % % % % % % % % % 
%         channel_h = [channel_h_coefficient zeros(1,length(signal_ori)-length(channel_h_coefficient))]; 
%         cloumn = [channel_h(1) zeros(1,length(y)-1)];
%         toep_h = toeplitz(cloumn,channel_h);
%         toep_h_inv = pinv(toep_h);
%         zf_filter = fliplr(toep_h_inv(1,:));
%         x_hat = conv(y,zf_filter);
%         % % % % % % % % % % % ZF end % % % % % % % % % % %
%         x_hat = x_hat(min(length(y),length(zf_filter)):end);
%         x_hat = x_hat./norm(x_hat,2)*sqrt(length(x_hat))*sqrt(5); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
% 
%         x_demod = (x_hat+3)/2;
%         x_demod(x_demod <= 0.5) = 0;
%         x_demod(0.5 < x_demod & x_demod <= 1.5) = 1;
%         x_demod(1.5 < x_demod & x_demod <= 2.5) = 2;
%         x_demod(2.5 < x_demod) = 3;
% 
%     elseif strcmpi(p.Results.demodulation,'MMSE')
%         r_phase = mod(index2,times);        % Find the phase of maximum correlation point
%         r_for_judge = r3(r_phase:times:end); 
%         [corr_degree,~] = sort(r_for_judge,'descend');
%         channel_h_coefficient = corr_degree(1:30);  % Find 30 points with maximum correlation as h1 ... h30
%         y = current_signal(index3:times:index3+times*(length(signal_ori)-1)); % Received signal represented by symbols.
% 
%         channel_h = [channel_h_coefficient zeros(1,length(signal_ori)-length(channel_h_coefficient))]; 
%         cloumn = [channel_h(1) zeros(1,length(y)-1)];
%         H = toeplitz(cloumn,channel_h);                                   % H
% 
%         G =  H'*inv((H*H'+sigma2/ps));
%         x_hat_t = G*y.';
%         x_hat = x_hat_t.';
%         x_hat = x_hat./norm(x_hat,2)*sqrt(length(x_hat))*sqrt(5);
% 
%         x_demod = (x_hat+3)/2;
%         x_demod(x_demod <= 0.5) = 0;
%         x_demod(0.5 < x_demod & x_demod <= 1.5) = 1;
%         x_demod(1.5 < x_demod & x_demod <= 2.5) = 2;
%         x_demod(2.5 < x_demod) = 3;
% 
%     end
% 
% end

