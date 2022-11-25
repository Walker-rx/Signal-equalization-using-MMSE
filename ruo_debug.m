clear;
close all;

channel_choice = 1; % % 1 is direct channel
channel_choice_inf = 4;
dir_up = "./data_set_final/";
test = [1 0 1 0 1 1 1];
% bias = input('bias(mA): ')

% dp821A = visa('ni','USB0::0x1AB1::0x0E11::DP8E163250125::INSTR');
% % dp821A = visa('ni','ASRL3:INSTR');
% fopen(dp821A);
% % fprintf(dp821A,'*IDN?'); 
% % aa = fscanf(dp821A);
% % display(aa);
% fprintf(dp821A,':OUTP:OCP:VAL CH2,1.500');
% fprintf(dp821A,':OUTP:OCP:VAL? CH2');
% aa = fscanf(dp821A);
% display(aa);
% fprintf(dp821A,':OUTP:OCP CH2,ON');
% fprintf(dp821A,':OUTP CH2,ON');
% fprintf(dp821A,':APPL CH1,24,1.000');
% fprintf(dp821A,':OUTP CH1,ON');

bias_name = 780;
pause(0.5);

M = 4;
data_length = 10000;
zero_length = 3000;
ls_order = 50;

origin_rate = 25e6; 
bw = origin_rate/2;         % baseband bandwidth
f_rate = 160e6;
d_rate = 150e6;

times = d_rate/origin_rate;
num_of_windows = 100;

% pilot = pilot_gen([0 0 0 1 1 1 1]);
% pilot = pilot_gen([1 1 1 0 0 0 1 1 1]);
pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 0 1]);    %  2^11 = 2048
% pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 1 0 1 1]);   %  2^13 = 8192
pilot_ps = bandpower(pilot);
pilot_length = length(pilot);
pilot_bpsk = pilot*2-1;

filter_ord = 1000;     % filter order used in function sam_rate_con
rp = 0.00057565;      
rst = 1e-4;       % filter parameter used in function sam_rate_con

llcm_transmit = lcm(origin_rate,f_rate);
upf_transmit = llcm_transmit/origin_rate;  % upsampling parameters
dof_transmit = llcm_transmit/f_rate;
ups_rate_transmit = origin_rate*upf_transmit;
filter_transmit = firceqrip(filter_ord,bw/(ups_rate_transmit/2),[rp rst],'passedge'); % Impulse function
filter_transmit = filter_transmit./norm(filter_transmit,2)*sqrt(bw/ups_rate_transmit*2); % filter used in function sam_rate_con in transmit part

llcm_receive = lcm(f_rate,d_rate);
upf_receive = llcm_receive/f_rate;  % upsampling parameters
dof_receive = llcm_receive/d_rate;
ups_rate_receive = f_rate*upf_receive;
filter_receive = firceqrip(filter_ord,bw/(ups_rate_receive/2),[rp rst],'passedge'); % Impulse function
filter_receive = filter_receive./norm(filter_receive,2)*sqrt(bw/ups_rate_receive*2); % filter used in function sam_rate_con in receive part

% origin_rate = gpuArray(double(origin_rate));
% f_rate = gpuArray(double(f_rate));
% d_rate = gpuArray(double(d_rate));
% num_of_windows = gpuArray(double(num_of_windows));
% times = gpuArray(double(times));
% pilot_length = gpuArray(double(pilot_length));
% delay = gpuArray(double(delay));
% filter_transmit = gpuArray(double(filter_transmit));
% filter_receive = gpuArray(double(filter_receive));
% upf_transmit = gpuArray(double(upf_transmit));
% dof_transmit = gpuArray(double(dof_transmit));
% upf_receive = gpuArray(double(upf_receive));
% dof_receive = gpuArray(double(dof_receive));
% h_channel = gpuArray(double(h_channel));
% h_channel_delay = gpuArray(double(h_channel_delay));

t = datetime('now');
save_path = "snr_ser/direct/"+t.Month+"."+t.Day+'/debug';

% amp_begin = -6;
amp_begin = 50;
amp_end = 53;
fprintf('add zero,ls order=%d,pilot length=%d .\n',ls_order,pilot_length);
for amp = 50:amp_end
    save_path_mat = save_path+"/amp"+amp+"/mat";
    save_path_txt = save_path+"/amp"+amp+"/txt";
    if(~exist(save_path_mat,'dir'))
        mkdir(char(save_path_mat));
    end
    if(~exist(save_path_txt,'dir'))
        mkdir(char(save_path_txt));
    end
    
    looptime = 0;
    ps_aftercorrect = 0;
    pn_aftercorrect = 0;
    ps_beforecorrect = 0;
    pn_beforecorrect = 0;
    errornum_zf = 0;
    errornum_mmse = 0;
    errornum_ls_aftercorrect = 0;
    errornum_ls_beforecorrect = 0;
    total_length_beforecorrect = 0;
    total_length_aftercorrect = 0;
    txerror_num = 0;
    fprintf('amp = %d .\n', amp);
    
    while(errornum_ls_aftercorrect <= 30 || looptime < 100)
%     while(errornum_zf <= 100 || errornum_mmse <= 100 || looptime < 50)
        
        looptime = looptime+1;
        
        data = randi([0,M-1],[1,data_length]);
        data_mpam = real(pammod(data,M));
        data_mpam_ps = bandpower(data_mpam);
        
        unmod = [pilot,zeros(1,zero_length),data];
        signal_ori = [pilot_bpsk,zeros(1,zero_length),data_mpam];
        signal_upsample = ruo_sam_rate_con(signal_ori,filter_transmit,upf_transmit,dof_transmit);       
        signal_send_tmp = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp);
        signal_send_tmp_inf = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp_end);
        signal_rand = rand([1,512]);
        signal_send = [signal_rand zeros(1,10) signal_send_tmp];
        signal_send_inf = [signal_rand zeros(1,10) signal_send_tmp_inf];
      
        ruo_pam4_send;
        ruo_pam4_send_correct;
        
        %% Locating the transmission error point
        signal_received_inf = ruo_sam_rate_con(signal_pass_channel_inf,filter_receive,upf_receive,dof_receive);
        
        [fin_syn_point_inf,coar_syn_point_inf] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_inf,num_of_windows);
        
        signal_demod_ls_inf = ruo_signal_equal(signal_ori,signal_received_inf,times,coar_syn_point_inf,pilot_length,zero_length,num_of_windows,{'ls',ls_order});
        
        data_demod_ls_inf = signal_demod_ls_inf(pilot_length+zero_length+1:end);
        
        if length(signal_demod_ls_inf) < length(unmod)
            compare_length_inf = length(signal_demod_ls_inf);
            error_location_inf = find(signal_demod_ls_inf ~= unmod(1:compare_length_inf));
        else
            compare_length_inf = length(unmod);
            error_location_inf = find(signal_demod_ls_inf(1:compare_length_inf) ~= unmod);
        end
        txerror_num_loop = length(error_location_inf) ;
        txerror_num = txerror_num + txerror_num_loop ;
        txerror_num_average = txerror_num/looptime ;
        %% Calculating the SER without correction
        signal_received_channel1_send1 = ruo_sam_rate_con(signal_pass_channel,filter_receive,upf_receive,dof_receive);       
        [fin_syn_point_channel1_send1 , coar_syn_point_channel1_send1] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_channel1_send1,num_of_windows);
        
        [length_loop_channel1_send1 , ps_loop_channel1_send1 , pn_loop_channel1_send1 , ...
                       errornum_ls_loop_channel1_send1 , error_location_loop_channel1_send1 , data_demod_ls_channel1_send1] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_channel1_send1 , pilot_length , zero_length , data_length ...
                                , coar_syn_point_channel1_send1 , fin_syn_point_channel1_send1 , times , num_of_windows , ls_order);
        
                            
        signal_received_channel1_send2 = ruo_sam_rate_con(signal_pass_channel_correct,filter_receive,upf_receive,dof_receive);
        [fin_syn_point_channel1_send2 , coar_syn_point_channel1_send2] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_channel1_send2,num_of_windows);
        
        [length_loop_channel1_send2 , ps_loop_channel1_send2 , pn_loop_channel1_send2 , ...
                      errornum_ls_loop_channel1_send2 , error_location_loop_channel1_send2 , data_demod_ls_channel1_send2] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_channel1_send2 , pilot_length , zero_length , data_length ...
                                , coar_syn_point_channel1_send2 , fin_syn_point_channel1_send2 , times , num_of_windows , ls_order);
        
        total_length_beforecorrect = total_length_beforecorrect + length_loop_channel1_send1 + length_loop_channel1_send2;
        ps_beforecorrect = ps_beforecorrect + ps_loop_channel1_send1 + ps_loop_channel1_send2;
        pn_beforecorrect = pn_beforecorrect + pn_loop_channel1_send1 + pn_loop_channel1_send2;
        errornum_ls_beforecorrect = errornum_ls_beforecorrect + errornum_ls_loop_channel1_send1 + errornum_ls_loop_channel1_send2;
        
        snr_ls_beforecorrect = 10*log10(ps_beforecorrect/pn_beforecorrect);
        ser_ls_beforecorrect = errornum_ls_beforecorrect/total_length_beforecorrect;
        
        %% Replacing the transmission error point
        signal_received_tmp = signal_received_channel1_send1;
        signal_received_aftercorrect = signal_received_channel1_send1;
        replace_loc = error_location_inf*times + fin_syn_point_channel1_send1 - 1;
        replace_loc_correct = error_location_inf*times + fin_syn_point_channel1_send2 - 1;
        replace_length = 8*times;
        for i = 1:txerror_num_loop
            signal_received_aftercorrect(replace_loc(i)-replace_length : replace_loc(i)+replace_length) = signal_received_channel1_send2(replace_loc_correct(i)-replace_length : replace_loc_correct(i)+replace_length);
        end       
        %% Calculating the SER after correction
        [length_loop_aftercorrect , ps_loop_aftercorrect , pn_loop_aftercorrect , ...
                    errornum_ls_loop_aftercorrect , error_location_aftercorrect , data_demod_ls_aftercorrect] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_aftercorrect , pilot_length , zero_length , data_length ...
                                , coar_syn_point_channel1_send1 , fin_syn_point_channel1_send1 , times , num_of_windows , ls_order);
        
        total_length_aftercorrect = total_length_aftercorrect + length_loop_aftercorrect;
        errornum_ls_aftercorrect = errornum_ls_aftercorrect + errornum_ls_loop_aftercorrect;
        ps_aftercorrect = ps_aftercorrect + ps_loop_aftercorrect;
        pn_aftercorrect = pn_aftercorrect + pn_loop_aftercorrect;
        
        ser_ls_aftercorrect = errornum_ls_aftercorrect/total_length_aftercorrect;
        snr_ls_aftercorrect = 10*log10(ps_aftercorrect/pn_aftercorrect);
        
%         snr_ls_beforecorrect = 0;
%         ser_ls_beforecorrect = 0; 
        if mod(looptime,5) == 0
           fprintf(' amp = %d , %d times, average tx error num = %d .\n',amp,looptime,txerror_num_average);
           fprintf(' %f times, before, snr = %d , total ls error num = %d , ls error rate = %.6g .\n',looptime,snr_ls_beforecorrect,errornum_ls_beforecorrect/2,ser_ls_beforecorrect);
           fprintf(' %f times, after, snr = %d , total ls error num = %d , ls error rate = %.6g .\n',looptime,snr_ls_aftercorrect,errornum_ls_aftercorrect,ser_ls_aftercorrect);
           fprintf(' %f times, before error = %d , replace num = %d , after error = %d .\n',looptime , errornum_ls_loop_channel1_send1 , txerror_num_loop ,  errornum_ls_loop_aftercorrect);
%            pause(5);
        end
    %% Save data     
        save_data = ['save_data_' num2str(looptime)];
        save_data_inf = ['save_data_inf_' num2str(looptime)];
        save_data_beforecorrect = ['save_data_beforecorrect_' num2str(looptime)];
        save_data_forcorrect = ['save_data_forcorrect_' num2str(looptime)];
        save_data_aftercorrect = ['save_data_aftercorrect_' num2str(looptime)]; 
        save_signal_send = ['save_signal_send_' num2str(looptime)];
        save_signal_ori = ['save_signal_ori_' num2str(looptime)];
        save_signal_received_inf = ['save_signal_received_inf_' num2str(looptime)];
        save_fin_syn_point_inf = ['save_fin_syn_point_inf_' num2str(looptime)];
        save_coar_syn_point_inf = ['save_coar_syn_point_inf_' num2str(looptime)];
        save_signal_received_channel1_send1 = ['save_signal_received_channel1_send1_' num2str(looptime)];
        save_fin_syn_point_channel1_send1 = ['save_fin_syn_point_channel1_send1_' num2str(looptime)];
        save_coar_syn_point_channel1_send1 = ['save_coar_syn_point_channel1_send1_' num2str(looptime)];
        save_signal_received_channel1_send2 = ['save_signal_received_channel1_send2_' num2str(looptime)];
        save_fin_syn_point_channel1_send2 = ['save_fin_syn_point_channel1_send2_' num2str(looptime)];
        save_coar_syn_point_channel1_send2 = ['save_coar_syn_point_channel1_send2_' num2str(looptime)];       
        save_signal_received_aftercorrect = ['save_signal_received_aftercorrect_' num2str(looptime)];
        save_replace_location = ['save_replace_location_' num2str(looptime)];
        save_errlocation_before = ['save_errlocation_before_' num2str(looptime)];
        save_errlocation_forcorrect = ['save_errlocation_forcorrect_' num2str(looptime)];
        save_errlocation_after = ['save_errlocation_after_' num2str(looptime)];
        
        eval([save_data,'=data;']);
        eval([save_data_inf,'=data_demod_ls_inf;']);
        eval([save_data_beforecorrect,'=data_demod_ls_channel1_send1;']);
        eval([save_data_forcorrect,'=data_demod_ls_channel1_send2;']);
        eval([save_data_aftercorrect,'=data_demod_ls_aftercorrect;']);       
        eval([save_signal_send,'=signal_send;']);
        eval([save_signal_ori,'=signal_ori;']);
        eval([save_signal_received_inf,'=signal_received_inf;']);
        eval([save_fin_syn_point_inf,'=fin_syn_point_inf;']);
        eval([save_coar_syn_point_inf,'=coar_syn_point_inf;']);
        eval([save_signal_received_channel1_send1,'=signal_received_channel1_send1;']);
        eval([save_fin_syn_point_channel1_send1,'=fin_syn_point_channel1_send1;']);
        eval([save_coar_syn_point_channel1_send1,'=coar_syn_point_channel1_send1;']);
        eval([save_signal_received_channel1_send2,'=signal_received_channel1_send2;']);
        eval([save_fin_syn_point_channel1_send2,'=fin_syn_point_channel1_send2;']);
        eval([save_coar_syn_point_channel1_send2,'=coar_syn_point_channel1_send2;']);        
        eval([save_signal_received_aftercorrect,'=signal_received_aftercorrect;']);
        eval([save_replace_location,'=error_location_inf;']);
        eval([save_errlocation_before,'=error_location_loop_channel1_send1;']);
        eval([save_errlocation_forcorrect,'=error_location_loop_channel1_send2;']);
        eval([save_errlocation_after,'=error_location_aftercorrect;']);
        
        if looptime == 1
            save(save_path_mat+"/save_data.mat",save_data);
            save(save_path_mat+"/save_data_inf.mat",save_data_inf);
            save(save_path_mat+"/save_data_beforecorrect.mat",save_data_beforecorrect);
            save(save_path_mat+"/save_data_forcorrect.mat",save_data_forcorrect);
            save(save_path_mat+"/save_data_aftercorrect.mat",save_data_aftercorrect);            
            save(save_path_mat+"/save_signal_send.mat",save_signal_send);
            save(save_path_mat+"/save_signal_ori.mat",save_signal_ori);
            save(save_path_mat+"/save_signal_received_inf.mat",save_signal_received_inf);
            save(save_path_mat+"/save_fin_syn_point_inf.mat",save_fin_syn_point_inf);
            save(save_path_mat+"/save_coar_syn_point_inf.mat",save_coar_syn_point_inf);
            save(save_path_mat+"/save_signal_received_channel1_send1.mat",save_signal_received_channel1_send1);
            save(save_path_mat+"/save_fin_syn_point_channel1_send1.mat",save_fin_syn_point_channel1_send1);
            save(save_path_mat+"/save_coar_syn_point_channel1_send1.mat",save_coar_syn_point_channel1_send1);
            save(save_path_mat+"/save_signal_received_channel1_send2.mat",save_signal_received_channel1_send2);
            save(save_path_mat+"/save_fin_syn_point_channel1_send2.mat",save_fin_syn_point_channel1_send2);
            save(save_path_mat+"/save_coar_syn_point_channel1_send2.mat",save_coar_syn_point_channel1_send2);
            save(save_path_mat+"/save_signal_received_aftercorrect.mat",save_signal_received_aftercorrect);
            save(save_path_mat+"/save_replace_location.mat",save_replace_location);
            save(save_path_mat+"/save_errlocation_before.mat",save_errlocation_before);
            save(save_path_mat+"/save_errlocation_forcorrect.mat",save_errlocation_forcorrect);
            save(save_path_mat+"/save_errlocation_after.mat",save_errlocation_after);
            
            save_errnum = fopen(save_path_txt+"/save_errnum.txt",'w');           
            fprintf(save_errnum,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
        else
            save(save_path_mat+"/save_data.mat",save_data,'-append');
            save(save_path_mat+"/save_data_inf.mat",save_data_inf,'-append');
            save(save_path_mat+"/save_data_beforecorrect.mat",save_data_beforecorrect,'-append');
            save(save_path_mat+"/save_data_forcorrect.mat",save_data_forcorrect,'-append');
            save(save_path_mat+"/save_data_aftercorrect.mat",save_data_aftercorrect,'-append');
            save(save_path_mat+"/save_signal_send.mat",save_signal_send,'-append');
            save(save_path_mat+"/save_signal_ori.mat",save_signal_ori,'-append');
            save(save_path_mat+"/save_signal_received_inf.mat",save_signal_received_inf,'-append');
            save(save_path_mat+"/save_fin_syn_point_inf.mat",save_fin_syn_point_inf,'-append');
            save(save_path_mat+"/save_coar_syn_point_inf.mat",save_coar_syn_point_inf,'-append');
            save(save_path_mat+"/save_signal_received_channel1_send1.mat",save_signal_received_channel1_send1,'-append');
            save(save_path_mat+"/save_fin_syn_point_channel1_send1.mat",save_fin_syn_point_channel1_send1,'-append');
            save(save_path_mat+"/save_coar_syn_point_channel1_send1.mat",save_coar_syn_point_channel1_send1,'-append');
            save(save_path_mat+"/save_signal_received_channel1_send2.mat",save_signal_received_channel1_send2,'-append');
            save(save_path_mat+"/save_fin_syn_point_channel1_send2.mat",save_fin_syn_point_channel1_send2,'-append');
            save(save_path_mat+"/save_coar_syn_point_channel1_send2.mat",save_coar_syn_point_channel1_send2,'-append');
            save(save_path_mat+"/save_signal_received_aftercorrect.mat",save_signal_received_aftercorrect,'-append');
            save(save_path_mat+"/save_replace_location.mat",save_replace_location,'-append');
            save(save_path_mat+"/save_errlocation_before.mat",save_errlocation_before,'-append');
            save(save_path_mat+"/save_errlocation_forcorrect.mat",save_errlocation_forcorrect,'-append');
            save(save_path_mat+"/save_errlocation_after.mat",save_errlocation_after,'-append');
            
            save_errnum = fopen(save_path_txt+"/save_errnum.txt",'a');
        end
        fprintf(save_errnum,' %f times, before error = %d , replace num = %d , after error = %d .\n', ...
                              looptime , errornum_ls_loop_channel1_send1 , txerror_num_loop ,  errornum_ls_loop_aftercorrect);
        fclose(save_errnum);
    end  
    %%
%     ser_ls = gather(ser_ls);

    %% Save data              
%     if amp == amp_begin
%         fsnr_beforecorrect = fopen(save_path+"/snr_beforecorrect.txt",'w');
%         fser_beforecorrect = fopen(save_path+"/ser_beforecorrect.txt",'w');
%         fsnr_aftercorrect = fopen(save_path+"/snr_aftercorrect.txt",'w');
%         fser_aftercorrect = fopen(save_path+"/ser_aftercorrect.txt",'w');
%         
%         fprintf(fsnr_beforecorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
%         fprintf(fser_beforecorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
%         fprintf(fsnr_aftercorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
%         fprintf(fser_aftercorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
%     else
%         fsnr_beforecorrect = fopen(save_path+"/snr_beforecorrect.txt",'a');
%         fser_beforecorrect = fopen(save_path+"/ser_beforecorrect.txt",'a');
%         fsnr_aftercorrect = fopen(save_path+"/snr_aftercorrect.txt",'a');
%         fser_aftercorrect = fopen(save_path+"/ser_aftercorrect.txt",'a');
%     end
%     fprintf(fsnr_beforecorrect,'%.8f \r\n',snr_ls_beforecorrect);
%     fprintf(fser_beforecorrect,'%.6g \r\n',ser_ls_beforecorrect);
%     fprintf(fsnr_aftercorrect,'%.8f \r\n',snr_ls_aftercorrect);
%     fprintf(fser_aftercorrect,'%.6g \r\n',ser_ls_aftercorrect);
%     fclose(fsnr_aftercorrect);
%     fclose(fser_aftercorrect);

end


