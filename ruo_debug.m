clear;
close all;

channel_choice = 1; %  1 is light channel
channel_choice_inf = 3;%  3 is direct channel
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
data_length_initial = 10000;
zero_length = 3000;
zero_length_forsyn = 200;
ls_order = 50;
num_of_windows = 100;

times = 6;
origin_rate_tmp = 10e6;
f_rate = 160e6;
d_rate_tmp = origin_rate_tmp*times;

filter_order = 1000;     % filter order used in function sam_rate_con
rp = 0.00057565;      
rst = 1e-4;       % filter parameter used in function sam_rate_con

%% Generate pilot
% pilot = pilot_gen([0 0 0 1 1 1 1]);
% pilot = pilot_gen([1 1 1 0 0 0 1 1 1]);
pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 0 1]);    %  2^11 = 2048
% pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 1 0 1 1]);   %  2^13 = 8192
pilot_ps = bandpower(pilot);
pilot_length = length(pilot);
pilot_bpsk = pilot*2-1;

%% Generate filter
[origin_rate , d_rate , upf_transmit , dof_transmit , filter_transmit , upf_receive , dof_receive , filter_receive] = ruo_filter_gen(origin_rate_tmp , f_rate , times , filter_order , rp , rst);
rate_change = 100*abs(origin_rate-origin_rate_tmp)/origin_rate_tmp;
fprintf("rate's change = %.8f%% \n",rate_change);
ups_time = upf_transmit/dof_transmit;
data_length = round(data_length_initial*8/ups_time);

%%
t = datetime('now');
save_path = "vol_save/"+t.Month+"."+t.Day+"/"+origin_rate_tmp/1e6+"M";
if(~exist(save_path,'dir'))
    mkdir(char(save_path));
end

amp_begin = 8;
amp_end = 62;
amp_inf = 38;

fprintf('add zero,ls order=%d,pilot length=%d .\n',ls_order,pilot_length);
for amp = 40:amp_end
    save_path_mat = save_path+"/amp"+amp+"/mat";
    save_path_txt = save_path+"/amp"+amp+"/txt";
    
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
    replace_length = 8*times;
    replace_valid_num = 0;
    replace_correct_num = 0;
    save_time = 0;

    while(errornum_ls_aftercorrect <= 100 || looptime < 500)
%     while(looptime < 10)
        
        looptime = looptime+1;
         %% Signal send
        data = randi([0,M-1],[1,data_length]);
        data_mpam = real(pammod(data,M));
        data_mpam_ps = bandpower(data_mpam);  
        
        unmod = [pilot,zeros(1,zero_length),data];
        signal_ori = [pilot_bpsk,zeros(1,zero_length),data_mpam];
        pilot_bpsk_forsyn = [pilot_bpsk,zeros(1,zero_length_forsyn)];
        
        pilot_upsample_forsyn = ruo_sam_rate_con(pilot_bpsk_forsyn,filter_transmit,upf_transmit,dof_transmit);
        signal_upsample = ruo_sam_rate_con(signal_ori,filter_transmit,upf_transmit,dof_transmit);

        pilot_send_forsyn = pilot_upsample_forsyn./norm(pilot_upsample_forsyn,2)*sqrt(length(pilot_upsample_forsyn))*100*1.1^(amp_inf);
        signal_send_tmp = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp);
        signal_send_tmp_inf = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp_inf);
       
        signal_rand = rand([1,512]);
        signal_send = [signal_rand pilot_send_forsyn signal_send_tmp];
        signal_send_inf = [signal_rand pilot_send_forsyn signal_send_tmp_inf];
      
        ruo_pam4_send;
        ruo_pam4_send_correct;
        
        %% Locating the transmission error point
        signal_received_inf = ruo_sam_rate_con(signal_pass_channel_inf,filter_receive,upf_receive,dof_receive);
        
        [fin_syn_point_inf_tmp,coar_syn_point_inf_tmp] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_inf,num_of_windows);
        
        if (length(signal_received_inf) - length(signal_received_inf(fin_syn_point_inf_tmp:end))) > (pilot_length + zero_length_forsyn)*times
            coar_syn_point_inf = coar_syn_point_inf_tmp;
            fin_syn_point_inf = fin_syn_point_inf_tmp;
        else
            coar_syn_point_inf = coar_syn_point_inf_tmp + pilot_length + zero_length_forsyn;
            fin_syn_point_inf = fin_syn_point_inf_tmp + (pilot_length + zero_length_forsyn)*times;
        end
               
        signal_demod_ls_inf = ruo_signal_equal_ls(signal_ori,signal_received_inf,times,fin_syn_point_inf,pilot_length,zero_length,ls_order);
  
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
        signal_received_real_send1 = ruo_sam_rate_con(signal_pass_channel,filter_receive,upf_receive,dof_receive);       
        [fin_syn_point_real_send1_tmp , coar_syn_point_real_send1_tmp] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_real_send1,num_of_windows);      
        if (length(signal_received_real_send1) - length(signal_received_real_send1(fin_syn_point_real_send1_tmp:end))) > (pilot_length + zero_length_forsyn)*times
            coar_syn_point_real_send1 = coar_syn_point_real_send1_tmp;
            fin_syn_point_real_send1 = fin_syn_point_real_send1_tmp;
        else
            coar_syn_point_real_send1 = coar_syn_point_real_send1_tmp + pilot_length + zero_length_forsyn;
            fin_syn_point_real_send1 = fin_syn_point_real_send1_tmp + (pilot_length + zero_length_forsyn)*times;
        end              
        [length_loop_real_send1 , ps_loop_real_send1 , pn_loop_real_send1 , ...
            errornum_ls_loop_real_send1 , error_location_loop_real_send1 , data_demod_ls_real_send1] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_real_send1 , pilot_length , zero_length , data_length ...
                                 , fin_syn_point_real_send1 , times , ls_order);
                                   
        signal_received_real_send2 = ruo_sam_rate_con(signal_pass_channel_correct,filter_receive,upf_receive,dof_receive);
        [fin_syn_point_real_send2_tmp , coar_syn_point_real_send2_tmp] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_real_send2,num_of_windows);
        if (length(signal_received_real_send2) - length(signal_received_real_send2(fin_syn_point_real_send2_tmp:end))) > (pilot_length + zero_length_forsyn)*times
            coar_syn_point_real_send2 = coar_syn_point_real_send2_tmp;
            fin_syn_point_real_send2 = fin_syn_point_real_send2_tmp;
        else
            coar_syn_point_real_send2 = coar_syn_point_real_send2_tmp + pilot_length + zero_length_forsyn;
            fin_syn_point_real_send2 = fin_syn_point_real_send2_tmp + (pilot_length + zero_length_forsyn)*times;
        end        
        [length_loop_real_send2 , ps_loop_real_send2 , pn_loop_real_send2 , ...
            errornum_ls_loop_real_send2 , error_location_loop_real_send2 , data_demod_ls_real_send2] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_real_send2 , pilot_length , zero_length , data_length ...
                                 , fin_syn_point_real_send2 , times , ls_order);
        
        total_length_beforecorrect = total_length_beforecorrect + length_loop_real_send1 + length_loop_real_send2;
        ps_beforecorrect = ps_beforecorrect + ps_loop_real_send1 + ps_loop_real_send2;
        pn_beforecorrect = pn_beforecorrect + pn_loop_real_send1 + pn_loop_real_send2;
        errornum_ls_beforecorrect = errornum_ls_beforecorrect + errornum_ls_loop_real_send1 + errornum_ls_loop_real_send2;
        
        snr_ls_beforecorrect = 10*log10(ps_beforecorrect/pn_beforecorrect);
        ser_ls_beforecorrect = errornum_ls_beforecorrect/total_length_beforecorrect;
        
        %% Replacing the transmission error point
        signal_received_tmp = signal_received_real_send1;
        signal_received_aftercorrect = signal_received_real_send1;
        signal_ori_resyn = signal_received_aftercorrect(fin_syn_point_real_send1_tmp:end);
        fin_syn_point_forcorrect_tmp = ruo_signal_syn_recorrect(signal_ori_resyn,signal_received_real_send2,fin_syn_point_real_send2_tmp);
        if fin_syn_point_real_send1 == fin_syn_point_real_send1_tmp && fin_syn_point_real_send2 == fin_syn_point_real_send2_tmp
            fin_syn_point_forcorrect = fin_syn_point_forcorrect_tmp;
        else         
            fin_syn_point_forcorrect = fin_syn_point_forcorrect_tmp + (pilot_length + zero_length_forsyn)*times;
        end
        if txerror_num_loop > 0
            replace_loc = error_location_inf*times + fin_syn_point_real_send1 - 1;
            replace_loc_correct = error_location_inf*times + fin_syn_point_forcorrect - 1;         
            for i = 1:txerror_num_loop
                if (replace_loc(i)-replace_length) > 0
                    if replace_loc(i)+replace_length > length(signal_received_aftercorrect) || replace_loc_correct(i)+replace_length > length(signal_received_real_send2)
                        signal_received_aftercorrect(length(signal_received_aftercorrect)-2*replace_length : end) = signal_received_real_send2(length(signal_received_real_send2)-2*replace_length : end);
                    else
                        signal_received_aftercorrect(replace_loc(i)-replace_length : replace_loc(i)+replace_length) = signal_received_real_send2(replace_loc_correct(i)-replace_length : replace_loc_correct(i)+replace_length);
                    end                   
                else
                    signal_received_aftercorrect(1 : replace_loc(i)+replace_length) = signal_received_real_send2(1 : replace_loc_correct(i)+replace_length);
                end
            end
        end
        %% Calculating the SER after correction
        [length_loop_aftercorrect , ps_loop_aftercorrect , pn_loop_aftercorrect , ...
            errornum_ls_loop_aftercorrect , error_location_aftercorrect , data_demod_ls_aftercorrect] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_aftercorrect , pilot_length , zero_length , data_length ...
                                 , fin_syn_point_real_send1 , times , ls_order);
        
        total_length_aftercorrect = total_length_aftercorrect + length_loop_aftercorrect;
        errornum_ls_aftercorrect = errornum_ls_aftercorrect + errornum_ls_loop_aftercorrect;
        ps_aftercorrect = ps_aftercorrect + ps_loop_aftercorrect;
        pn_aftercorrect = pn_aftercorrect + pn_loop_aftercorrect;
        
        ser_ls_aftercorrect = errornum_ls_aftercorrect/total_length_aftercorrect;
        snr_ls_aftercorrect = 10*log10(ps_aftercorrect/pn_aftercorrect);
        
        replace_correct_num_tmp = errornum_ls_loop_real_send1 - errornum_ls_loop_aftercorrect;
        replace_correct_num = replace_correct_num + replace_correct_num_tmp;
        if errornum_ls_loop_real_send1 > errornum_ls_loop_aftercorrect
            replace_valid_num = replace_valid_num + 1;
        end
%         snr_ls_beforecorrect = 0;
%         ser_ls_beforecorrect = 0; 
        fprintf(' %f times, before error = %d , replace num = %d , after error = %d , replace correct num = %d .\n',looptime , errornum_ls_loop_real_send1 , txerror_num_loop ,  errornum_ls_loop_aftercorrect , replace_correct_num_tmp);
        if mod(looptime,5) == 0
           fprintf('\n');
           fprintf(' amp = %d , %d times, average tx error num = %d .\n',amp,looptime,txerror_num_average);
           fprintf(' %f times, before, snr = %d , total ls error num = %d , ls error rate = %.6g .\n',looptime,snr_ls_beforecorrect,errornum_ls_beforecorrect/2,ser_ls_beforecorrect);
           fprintf(' %f times, after, snr = %d , total ls error num = %d , ls error rate = %.6g .\n',looptime,snr_ls_aftercorrect,errornum_ls_aftercorrect,ser_ls_aftercorrect);
           fprintf(' %f times, total replace correct num = %d .\n\n',looptime , replace_correct_num);
%            pause(1)
        end
                 
    %% Save data  
        if(~exist(save_path_mat,'dir'))
            mkdir(char(save_path_mat));
        end
        if(~exist(save_path_txt,'dir'))
            mkdir(char(save_path_txt));
        end

        save_time = save_time + 1;
        save_data = ['save_data_' num2str(looptime)];
        save_data_inf = ['save_data_inf_' num2str(looptime)];
        save_data_beforecorrect = ['save_data_beforecorrect_' num2str(looptime)];
        save_data_forcorrect = ['save_data_forcorrect_' num2str(looptime)];
        save_data_aftercorrect = ['save_data_aftercorrect_' num2str(looptime)];
        save_signal_send = ['save_signal_send_' num2str(looptime)];
        save_signal_ori = ['save_signal_ori_' num2str(looptime)];           
        save_signal_send_inf = ['save_signal_send_inf_' num2str(looptime)];
        save_signal_pass_channel_inf = ['save_signal_pass_channel_inf_' num2str(looptime)];
        save_signal_pass_channel = ['save_signal_pass_channel_' num2str(looptime)];
        save_signal_pass_channel_correct = ['save_signal_pass_channel_correct_' num2str(looptime)];
        save_signal_received_inf = ['save_signal_received_inf_' num2str(looptime)];
        save_fin_syn_point_inf = ['save_fin_syn_point_inf_' num2str(looptime)];
        save_coar_syn_point_inf = ['save_coar_syn_point_inf_' num2str(looptime)];
        save_signal_received_real_send1 = ['save_signal_received_real_send1_' num2str(looptime)];
        save_fin_syn_point_real_send1 = ['save_fin_syn_point_real_send1_' num2str(looptime)];
        save_coar_syn_point_real_send1 = ['save_coar_syn_point_real_send1_' num2str(looptime)];
        save_signal_received_real_send2 = ['save_signal_received_real_send2_' num2str(looptime)];
        save_fin_syn_point_forcorrect = ['save_fin_syn_point_forcorrect_' num2str(looptime)];
        save_coar_syn_point_real_send2 = ['save_coar_syn_point_real_send2_' num2str(looptime)];
        save_signal_received_aftercorrect = ['save_signal_received_aftercorrect_' num2str(looptime)];
        save_replace_location = ['save_replace_location_' num2str(looptime)];
        save_errlocation_before = ['save_errlocation_before_' num2str(looptime)];
        save_errlocation_forcorrect = ['save_errlocation_forcorrect_' num2str(looptime)];
        save_errlocation_after = ['save_errlocation_after_' num2str(looptime)];

        eval([save_data,'=data;']);
        eval([save_data_inf,'=data_demod_ls_inf;']);
        eval([save_data_beforecorrect,'=data_demod_ls_real_send1;']);
        eval([save_data_forcorrect,'=data_demod_ls_real_send2;']);
        eval([save_data_aftercorrect,'=data_demod_ls_aftercorrect;']);
        eval([save_signal_send,'=signal_send;']);
        eval([save_signal_ori,'=signal_ori;']);       
        eval([save_signal_send_inf,'=signal_send_inf;']);
        eval([save_signal_pass_channel_inf,'=signal_pass_channel_inf;']);
        eval([save_signal_pass_channel,'=signal_pass_channel;']);
        eval([save_signal_pass_channel_correct,'=signal_pass_channel_correct;']);       
        eval([save_signal_received_inf,'=signal_received_inf;']);
        eval([save_fin_syn_point_inf,'=fin_syn_point_inf;']);
        eval([save_coar_syn_point_inf,'=coar_syn_point_inf;']);
        eval([save_signal_received_real_send1,'=signal_received_real_send1;']);
        eval([save_fin_syn_point_real_send1,'=fin_syn_point_real_send1;']);
        eval([save_coar_syn_point_real_send1,'=coar_syn_point_real_send1;']);
        eval([save_signal_received_real_send2,'=signal_received_real_send2;']);
        eval([save_fin_syn_point_forcorrect,'=fin_syn_point_forcorrect;']);
        eval([save_coar_syn_point_real_send2,'=coar_syn_point_real_send2;']);
        eval([save_signal_received_aftercorrect,'=signal_received_aftercorrect;']);
        eval([save_replace_location,'=error_location_inf;']);
        eval([save_errlocation_before,'=error_location_loop_real_send1;']);
        eval([save_errlocation_forcorrect,'=error_location_loop_real_send2;']);
        eval([save_errlocation_after,'=error_location_aftercorrect;']);

        if save_time == 1
            save(save_path_mat+"/save_data.mat",save_data);
            save(save_path_mat+"/save_data_inf.mat",save_data_inf);
            save(save_path_mat+"/save_data_beforecorrect.mat",save_data_beforecorrect);
            save(save_path_mat+"/save_data_forcorrect.mat",save_data_forcorrect);
            save(save_path_mat+"/save_data_aftercorrect.mat",save_data_aftercorrect);
            save(save_path_mat+"/save_signal_send.mat",save_signal_send);
            save(save_path_mat+"/save_signal_ori.mat",save_signal_ori);            
            save(save_path_mat+"/save_signal_send_inf.mat",save_signal_send_inf);
            save(save_path_mat+"/save_signal_pass_channel_inf.mat",save_signal_pass_channel_inf);
            save(save_path_mat+"/save_signal_pass_channel.mat",save_signal_pass_channel);
            save(save_path_mat+"/save_signal_pass_channel_correct.mat",save_signal_pass_channel_correct);          
            save(save_path_mat+"/save_signal_received_inf.mat",save_signal_received_inf);
            save(save_path_mat+"/save_fin_syn_point_inf.mat",save_fin_syn_point_inf);
            save(save_path_mat+"/save_coar_syn_point_inf.mat",save_coar_syn_point_inf);
            save(save_path_mat+"/save_signal_received_real_send1.mat",save_signal_received_real_send1);
            save(save_path_mat+"/save_fin_syn_point_real_send1.mat",save_fin_syn_point_real_send1);
            save(save_path_mat+"/save_coar_syn_point_real_send1.mat",save_coar_syn_point_real_send1);
            save(save_path_mat+"/save_signal_received_real_send2.mat",save_signal_received_real_send2);
            save(save_path_mat+"/save_fin_syn_point_forcorrect.mat",save_fin_syn_point_forcorrect);
            save(save_path_mat+"/save_coar_syn_point_real_send2.mat",save_coar_syn_point_real_send2);
            save(save_path_mat+"/save_signal_received_aftercorrect.mat",save_signal_received_aftercorrect);
            save(save_path_mat+"/save_replace_location.mat",save_replace_location);
            save(save_path_mat+"/save_errlocation_before.mat",save_errlocation_before);
            save(save_path_mat+"/save_errlocation_forcorrect.mat",save_errlocation_forcorrect);
            save(save_path_mat+"/save_errlocation_after.mat",save_errlocation_after);
            
            save_errnum = fopen(save_path_txt+"/save_errnum.txt",'w');
            fprintf(save_errnum,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , origin_rate = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
        else
            save(save_path_mat+"/save_data.mat",save_data,'-append');
            save(save_path_mat+"/save_data_inf.mat",save_data_inf,'-append');
            save(save_path_mat+"/save_data_beforecorrect.mat",save_data_beforecorrect,'-append');
            save(save_path_mat+"/save_data_forcorrect.mat",save_data_forcorrect,'-append');
            save(save_path_mat+"/save_data_aftercorrect.mat",save_data_aftercorrect,'-append');
            save(save_path_mat+"/save_signal_send.mat",save_signal_send,'-append');
            save(save_path_mat+"/save_signal_ori.mat",save_signal_ori,'-append');           
            save(save_path_mat+"/save_signal_send_inf.mat",save_signal_send_inf,'-append');
            save(save_path_mat+"/save_signal_pass_channel_inf.mat",save_signal_pass_channel_inf,'-append');
            save(save_path_mat+"/save_signal_pass_channel.mat",save_signal_pass_channel,'-append');
            save(save_path_mat+"/save_signal_pass_channel_correct.mat",save_signal_pass_channel_correct,'-append');          
            save(save_path_mat+"/save_signal_received_inf.mat",save_signal_received_inf,'-append');
            save(save_path_mat+"/save_fin_syn_point_inf.mat",save_fin_syn_point_inf,'-append');
            save(save_path_mat+"/save_coar_syn_point_inf.mat",save_coar_syn_point_inf,'-append');
            save(save_path_mat+"/save_signal_received_real_send1.mat",save_signal_received_real_send1,'-append');
            save(save_path_mat+"/save_fin_syn_point_real_send1.mat",save_fin_syn_point_real_send1,'-append');
            save(save_path_mat+"/save_coar_syn_point_real_send1.mat",save_coar_syn_point_real_send1,'-append');
            save(save_path_mat+"/save_signal_received_real_send2.mat",save_signal_received_real_send2,'-append');
            save(save_path_mat+"/save_fin_syn_point_forcorrect.mat",save_fin_syn_point_forcorrect,'-append');
            save(save_path_mat+"/save_coar_syn_point_real_send2.mat",save_coar_syn_point_real_send2,'-append');
            save(save_path_mat+"/save_signal_received_aftercorrect.mat",save_signal_received_aftercorrect,'-append');
            save(save_path_mat+"/save_replace_location.mat",save_replace_location,'-append');
            save(save_path_mat+"/save_errlocation_before.mat",save_errlocation_before,'-append');
            save(save_path_mat+"/save_errlocation_forcorrect.mat",save_errlocation_forcorrect,'-append');
            save(save_path_mat+"/save_errlocation_after.mat",save_errlocation_after,'-append');

            save_errnum = fopen(save_path_txt+"/save_errnum.txt",'a');
        end
        fprintf(save_errnum,' %f times, before error = %d , replace num = %d , after error = %d , replace correct num = %d .\n', ...
            looptime , errornum_ls_loop_real_send1 , txerror_num_loop ,  errornum_ls_loop_aftercorrect , replace_correct_num_tmp);
        fclose(save_errnum);
    end
end