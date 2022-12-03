clear;
close all;

t = datetime('now');

channel_choice = 1;  % 1 is direct channel
channel_choice2 = 4;
dir_up = "./data_set_final/";
test = [1 0 1 0 1 1 1];

pilot_length = 2047;
data_length = 40000;
zero_length = 10000;
zero_length1 = 5000;
ls_order = 50;

origin_rate = 25e6; 
bw = origin_rate/2;         % baseband bandwidth
f_rate = 160e6;
d_rate = 150e6;

times = d_rate/origin_rate;
num_of_windows = 100;

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

amp_begin = 40;
amp_end = 52;
fprintf('add zero,ls order=%d,pilot length=%d .\n',ls_order,pilot_length);
for amp = amp_begin:amp_end
    save_path_voltest = "vol_save/11.12_test6/amp"+amp;
   if(~exist(save_path_voltest,'dir'))
        mkdir(char(save_path_voltest));
    end
    load("vol_save/11.4/amp"+amp+"/signal_ori_save.mat");
    looptime = 0;
    ps = 0;
    pn = 0;
    snr_sum = 0;
    errornum_ls = 0; 
    total_length = 0;
    fprintf('amp = %d .\n', amp);
    
%     while(errornum_ls <= 30 || looptime < 2000)
    while (looptime<=4)
        looptime = looptime+1;
        mat_location = looptime * 80;
        namestr_ori = ['signal_ori_save_' num2str(mat_location)];
        signal_ori = eval(namestr_ori);
        namestr_orisave = ['signal_save_ori_' num2str(mat_location)];
        eval([namestr_orisave,'= signal_ori;']);
        save_path = save_path_voltest+"/mat_location"+mat_location;
        
        pilot_bpsk = signal_ori(1:pilot_length);
        data_mpam = signal_ori(pilot_length+zero_length+1:pilot_length+zero_length+1+data_length-1);
        data = (data_mpam+3)/2;
        
        signal_ori2 = [pilot_bpsk,zeros(1,zero_length1),data_mpam];
        if(~exist(save_path,'dir'))
            mkdir(char(save_path));
        end
%         save(save_path+"/signal_save_ori"+"_amp"+amp+"_loc"+mat_location+".mat",namestr_orisave);
        exp_time = 0;
        false_num = 0;
        while (exp_time <=9)
            exp_time = exp_time + 1;
            
            ruo_pam4_testsend2;
            ruo_pam4_testsend;
            
            signal_received = ruo_sam_rate_con(signal_pass_channel,filter_receive,upf_receive,dof_receive);
            
            [fin_syn_point,coar_syn_point] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received,num_of_windows);
            
            signal_downsample = signal_received(fin_syn_point:times:end);         
            signal_fin = signal_received(fin_syn_point:end);
                     
            signal_demod_ls = ruo_signal_equal(signal_ori,signal_received,times,coar_syn_point,pilot_length,zero_length1,num_of_windows,{'ls',ls_order});
                      
            data_demod_ls = signal_demod_ls(pilot_length+zero_length1+1:end);
            
            if length(data_demod_ls) < data_length
                compare_length = length(data_demod_ls);
                error_location = find(data_demod_ls ~= data(1:compare_length));
            else
                compare_length = length(data);
                error_location = find(data_demod_ls(1:compare_length) ~= data);
            end           
            %%
            signal_received1 = ruo_sam_rate_con(signal_pass_channel1,filter_receive,upf_receive,dof_receive);
            
            [fin_syn_point1,coar_syn_point1] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received1,num_of_windows);
            
            signal_downsample1 = signal_received1(fin_syn_point1:times:end);        
            signal_fin1 = signal_received1(fin_syn_point1:end);
            
            signal_demod_ls1 = ruo_signal_equal(signal_ori,signal_received1,times,coar_syn_point1,pilot_length,zero_length1,num_of_windows,{'ls',ls_order});
            
            data_demod_ls1 = signal_demod_ls1(pilot_length+zero_length1+1:end);
            
            if length(data_demod_ls1) < data_length
                compare_length1 = length(data_demod_ls1);
                error_location1 = find(data_demod_ls1 ~= data(1:compare_length1));
            else
                compare_length1 = length(data);
                error_location1 = find(data_demod_ls1(1:compare_length1) ~= data);
            end
            
            %%
            signal_received2 = ruo_sam_rate_con(signal_pass_channel2,filter_receive,upf_receive,dof_receive);
            
            [fin_syn_point2,coar_syn_point2] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received2,num_of_windows);
            
            signal_downsample2 = signal_received2(fin_syn_point2:times:end);                     
            signal_fin2 = signal_received2(fin_syn_point2:end);
                      
            signal_demod_ls2 = ruo_signal_equal(signal_ori,signal_received2,times,coar_syn_point2,pilot_length,zero_length1,num_of_windows,{'ls',ls_order});
            
            data_demod_ls2 = signal_demod_ls2(pilot_length+zero_length1+1:end);
            
            if length(data_demod_ls2) < data_length
                compare_length2 = length(data_demod_ls2);
                error_location2 = find(data_demod_ls2 ~= data(1:compare_length2));
            else
                compare_length2 = length(data);
                error_location2 = find(data_demod_ls2(1:compare_length2) ~= data);
            end            
            %%
            inf_errloc = round((error_location2+zero_length1+pilot_length)*160/25);
            exp_errloc = round((error_location+zero_length1+pilot_length)*160/25);
            exp_errloc1 = round((error_location1+zero_length1+pilot_length)*160/25);
            
            [corr_point,corr_index,~] = intersect(exp_errloc,exp_errloc1);
            corr_num = length(corr_point);
            
            namestr_loc = ['errloc_save_' num2str(exp_time)];
            eval([namestr_loc,'= exp_errloc;']);
            namestr_loc1 = ['errloc_save1_' num2str(exp_time)];
            eval([namestr_loc1,'= exp_errloc1;']);
            namestr_corrindex = ['corrindex_save_' num2str(exp_time)];
            eval([namestr_corrindex,'= corr_index;']);
            
            inf_length = length(inf_errloc);
            exp_length = length(exp_errloc);
            exp_length1 = length(exp_errloc1);
            for i=1:length(inf_errloc)
                exp_errloc(find(inf_errloc(i)-20<=exp_errloc & exp_errloc<=inf_errloc(i)+20)) = [];
                exp_errloc1(find(inf_errloc(i)-20<=exp_errloc1 & exp_errloc1<=inf_errloc(i)+20)) = [];
            end
            undel_num = length(exp_errloc);
            undel_num1 = length(exp_errloc1);

            if corr_num ~=0
               false_num = false_num + 1; 
            end
%             namestr_loc = ['errloc_save_' num2str(exp_time)];
%             eval([namestr_loc,'= exp_errloc;']);
%             namestr_loc1 = ['errloc_save1_' num2str(exp_time)];
%             eval([namestr_loc1,'= exp_errloc1;']);
%             namestr_corrindex = ['corrindex_save_' num2str(exp_time)];
%             eval([namestr_corrindex,'= corr_index;']);
            if mod(exp_time,1) == 0
                fprintf(' amp = %d , %f times, exp_time = %d, .\n',amp,looptime,exp_time);
            end
            if exp_time == 1
                undelnum_save = fopen(save_path+"/unquit_num"+"_amp"+amp+"_loc"+mat_location+".txt",'w');
                save(save_path+"/errloc_save"+"_amp"+amp+"_loc"+mat_location+".mat",namestr_loc);
                save(save_path+"/errloc_save1"+"_amp"+amp+"_loc"+mat_location+".mat",namestr_loc1);
                save(save_path+"/corrindex_save"+"_amp"+amp+"_loc"+mat_location+".mat",namestr_corrindex);
            else
                undelnum_save = fopen(save_path+"/unquit_num"+"_amp"+amp+"_loc"+mat_location+".txt",'a');
                save(save_path+"/errloc_save"+"_amp"+amp+"_loc"+mat_location+".mat",namestr_loc,'-append');
                save(save_path+"/errloc_save1"+"_amp"+amp+"_loc"+mat_location+".mat",namestr_loc1,'-append');
                save(save_path+"/corrindex_save"+"_amp"+amp+"_loc"+mat_location+".mat",namestr_corrindex,'-append');
            end
            fprintf(undelnum_save,'send length = %d, another send length = %d, pilot = %d, zero = %d, data = %d \r\n',send_length,send_length2,pilot_length,zero_length1,data_length);
            fprintf(undelnum_save,'inf num = %d, exp num = %d, another exp num = %d, undelete num =%d, another undelete num =%d, corr num =%d, \r\n',inf_length,exp_length,exp_length1,undel_num,undel_num1,corr_num);
            fprintf(undelnum_save,'\r\n');
            if exp_time == 20
                fprintf(undelnum_save,'false time =%d, \r\n',false_num);
            end
            fclose(undelnum_save); 
            fprintf(' amp = %d , %f times, exp_time = %d, undelete num = %d, another undelete num = %d, corr num =%d .\n',amp,looptime,exp_time,undel_num,undel_num1,corr_num);
        end        
    end
end