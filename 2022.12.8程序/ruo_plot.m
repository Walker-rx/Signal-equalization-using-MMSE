clearvars -except filter_receive upf_receive dof_receive origin_rate d_rate num_of_windows times pilot_length zero_length num_of_windows ls_order

amp = 60;
looptime = 11;

ruo_load_data

signal_nofalse_tmp = ruo_sam_rate_con(signal_send,filter_receive,upf_receive,dof_receive);
[fin_syn_point_nofalse,coar_syn_point_nofalse] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_nofalse_tmp,num_of_windows);
signal_nofalse = signal_nofalse_tmp(fin_syn_point_nofalse:end);
% signal_nofalse = signal_nofalse/max(signal_nofalse);
signal_nofalse = signal_nofalse./norm(signal_nofalse,2)*sqrt(length(signal_nofalse));

signal_inf = signal_received_inf(fin_syn_point_inf:end)/max(signal_received_inf(fin_syn_point_inf:end));
% signal_notcorrect = signal_received_notcorrect(fin_syn_point:end)/max(signal_received_notcorrect(fin_syn_point:end));
signal_notcorrect = signal_received_notcorrect(fin_syn_point:end)./norm(signal_received_notcorrect(fin_syn_point:end),2)*sqrt(length(signal_received_notcorrect(fin_syn_point:end)));
signal_forcorrect = signal_received_forcorrect(fin_syn_point_forcorrect:end)/max(signal_received_forcorrect(fin_syn_point_forcorrect:end));
signal_aftercorrect = signal_received_aftercorrect(fin_syn_point:end)/max(signal_received_aftercorrect(fin_syn_point:end));

plot_point_replace = replace_location*times - 1;
plot_point_before = (errlocation_before+zero_length+pilot_length)*times - 1;
plot_point_forcorrect = (errlocation_forcorrect+zero_length+pilot_length)*times - 1;
plot_point_after = (errlocation_after+zero_length+pilot_length)*times - 1;
plot_point_data_before = errlocation_before - 1;
plot_point_data_forcorrect = errlocation_forcorrect - 1;
plot_point_data_after = errlocation_after - 1;

% %% �滻λ�õ�����
% for i = 1:length(plot_point_replace)
%     close all    
%     %%
%     figure
%     plot(signal_notcorrect,'b')
%     hold on
%     plot(signal_aftercorrect,'r+-')
%     plot(signal_forcorrect,'g')
%     axis(([plot_point_replace(i)-50 plot_point_replace(i)+50 -1.2 1.2]))
%     text(plot_point_replace(i),signal_notcorrect(plot_point_replace(i)),'o');
%     title("�滻��λ��");
%     legend('δ����������','�����������','��������������');    
%     %%    
%     figure
%     plot(signal_notcorrect,'b')
%     hold on
%     plot(signal_aftercorrect,'r')
%     axis(([plot_point_replace(i)-50 plot_point_replace(i)+50 -1.2 1.2]))
%     text(plot_point_replace(i),signal_notcorrect(plot_point_replace(i)),'o');
%     title("�滻��λ��");
%     legend('δ����������','�����������');   
%     %%
%     figure
%     plot(signal_nofalse,'b')
%     hold on
%     plot(signal_notcorrect,'r')
%     plot(signal_aftercorrect,'g')
%     axis(([plot_point_replace(i)-50 plot_point_replace(i)+50 -1.2 1.2]))
%     text(plot_point_replace(i),signal_notcorrect(plot_point_replace(i)),'o');
%     title("�滻��λ��");
%     legend('δ�����ŵ�������','δ����������','�����������');    
%     %%    
%     figure
%     plot(signal_nofalse,'b')
%     hold on
%     plot(signal_forcorrect,'r')
%     axis(([plot_point_replace(i)-50 plot_point_replace(i)+50 -1.2 1.2]))
%     text(plot_point_replace(i),signal_forcorrect(plot_point_replace(i)),'o');
%     title("�滻��λ��");
%     legend('δ�����ŵ�������','��������������');    
%     %%
%     
%     figure
%     plot(signal_nofalse,'b')
%     hold on
%     plot(signal_notcorrect,'r')
%     plot(signal_inf,'g')
%     axis(([plot_point_replace(i)-50 plot_point_replace(i)+50 -1.2 1.2]))
%     text(plot_point_replace(i),signal_notcorrect(plot_point_replace(i)),'o');
%     title("�滻��λ��");
%     legend('δ�����ŵ�������','δ����������','�ο��ŵ�����');  
% end
% 

% %% �滻ǰ�Ĵ���λ�õ�����
% for i = 1:length(plot_point_before)
%     close all
%     %% data
%     figure
%     plot(data_demod_beforecorrect,'b')
%     hold on
%     plot(data_demod_aftercorrect,'r+-')
%     plot(data_demod_forcorrect,'g')
%     axis(([plot_point_data_before(i)-50 plot_point_data_before(i)+50 -0.5 3.5]))
%     text(plot_point_data_before(i),signal_notcorrect(plot_point_data_before(i)),'o');
%     title("�滻��Ĵ���λ�õĽ��������");
%     legend('δ����������','�����������','��������������');
%     %%
%     figure
%     plot(data_demod_beforecorrect,'r')
%     hold on
%     plot(data_demod_inf,'g')
%     plot(data_demod_nofalse,'b+-')
%     axis(([plot_point_data_before(i)-50 plot_point_data_before(i)+50 -0.5 3.5]))
%     text(plot_point_data_before(i),signal_notcorrect(plot_point_data_before(i)),'o');
%     title("�滻��Ĵ���λ�õĽ��������");
%     legend('δ����������','�ο��ŵ�����','û�д������');
%     %%
%     figure
%     plot(data_demod_beforecorrect,'b')
%     hold on
%     plot(data_demod_nofalse,'r+-')
%     plot(data_demod_forcorrect,'g')
%     axis(([plot_point_data_before(i)-50 plot_point_data_before(i)+50 -0.5 3.5]))
%     text(plot_point_data_before(i),signal_notcorrect(plot_point_data_before(i)),'o');
%     title("�滻��Ĵ���λ�õĽ��������");
%     legend('δ����������','û�д������','��������������');
%     
%     %% signal
%     figure
%     plot(signal_notcorrect,'b')
%     hold on
%     plot(signal_aftercorrect,'r+-')
%     plot(signal_forcorrect,'g')
%     axis(([plot_point_before(i)-50 plot_point_before(i)+50 -1.2 1.2]))
%     text(plot_point_before(i),signal_notcorrect(plot_point_before(i)),'o');
%     title("�滻ǰ�Ĵ���λ��");
%     legend('δ����������','�����������','��������������');
%     %%
%     figure
%     plot(signal_notcorrect,'b')
%     hold on
%     plot(signal_aftercorrect,'r')
%     axis(([plot_point_before(i)-50 plot_point_before(i)+50 -1.2 1.2]))
%     text(plot_point_before(i),signal_notcorrect(plot_point_before(i)),'o');
%     title("�滻ǰ�Ĵ���λ��");
%     legend('δ����������','�����������');
%     %%
%     figure
%     plot(signal_nofalse,'b')
%     hold on
%     plot(signal_notcorrect,'r')
%     plot(signal_aftercorrect,'g')
%     axis(([plot_point_before(i)-50 plot_point_before(i)+50 -1.2 1.2]))
%     text(plot_point_before(i),signal_notcorrect(plot_point_before(i)),'o');
%     title("�滻ǰ�Ĵ���λ��");
%     legend('δ�����ŵ�������','δ����������','�����������');
%     %%
%     figure
%     plot(signal_nofalse,'b')
%     hold on
%     plot(signal_forcorrect,'r')
%     axis(([plot_point_before(i)-50 plot_point_before(i)+50 -1.2 1.2]))
%     text(plot_point_before(i),signal_forcorrect(plot_point_before(i)),'o');
%     title("�滻ǰ�Ĵ���λ��");
%     legend('δ�����ŵ�������','��������������');
%     %%
%     figure
%     plot(signal_nofalse,'b')
%     hold on
%     plot(signal_notcorrect,'r')
%     plot(signal_inf,'g')
%     axis(([plot_point_before(i)-50 plot_point_before(i)+50 -1.2 1.2]))
%     text(plot_point_before(i),signal_notcorrect(plot_point_before(i)),'o');
%     title("�滻ǰ�Ĵ���λ��");
%     legend('δ�����ŵ�������','δ����������','�ο��ŵ�����');
% end


%% �滻��Ĵ���λ�õ�����
for i = 200:length(plot_point_after)
    close all    
    %% data
    figure
    plot(data_demod_notcorrect,'b')
    hold on
    plot(data_demod_aftercorrect,'r+-')
    plot(data_demod_forcorrect,'g')
    axis(([plot_point_data_after(i)-50 plot_point_data_after(i)+50 -0.5 3.5]))
    text(plot_point_data_after(i),data_demod_aftercorrect(plot_point_data_after(i)),'o');
    title("�滻��Ĵ���λ�õĽ��������");
    legend('δ����������','�����������','��������������');
    %%
    figure
    plot(data_demod_notcorrect,'b')
    hold on
    plot(data_demod_aftercorrect,'r+-')
    axis(([plot_point_data_after(i)-50 plot_point_data_after(i)+50 -0.5 3.5]))
    text(plot_point_data_after(i),data_demod_aftercorrect(plot_point_data_after(i)),'o');
    title("�滻��Ĵ���λ�õĽ��������");
    legend('δ����������','�����������');
    %% 
    figure
    plot(data_demod_aftercorrect,'b')
    hold on
    plot(data_demod_forcorrect,'g')
    axis(([plot_point_data_after(i)-50 plot_point_data_after(i)+50 -0.5 3.5]))
    text(plot_point_data_after(i),data_demod_aftercorrect(plot_point_data_after(i)),'o');
    title("�滻��Ĵ���λ�õĽ��������");
    legend('�����������','��������������');
    %% 
    figure
    plot(data_demod_notcorrect,'b')
    hold on
    plot(data_demod_forcorrect,'g')
    axis(([plot_point_data_after(i)-50 plot_point_data_after(i)+50 -0.5 3.5]))
    text(plot_point_data_after(i),data_demod_forcorrect(plot_point_data_after(i)),'o');
    title("�滻��Ĵ���λ�õĽ��������");
    legend('δ����������','��������������');
    %%
    figure
    plot(data_demod_notcorrect,'r')
    hold on
    plot(data_demod_inf,'g')
    plot(data_demod_nofalse,'b+-')
    axis(([plot_point_data_after(i)-50 plot_point_data_after(i)+50 -0.5 3.5]))
    text(plot_point_data_after(i),data_demod_nofalse(plot_point_data_after(i)),'o');
    title("�滻��Ĵ���λ�õĽ��������");
    legend('δ����������','�ο��ŵ�����','û�д������');
    %%
    figure
    plot(data_demod_notcorrect,'b')
    hold on
    plot(data_demod_nofalse,'r+-')
    plot(data_demod_forcorrect,'g')
    axis(([plot_point_data_after(i)-50 plot_point_data_after(i)+50 -0.5 3.5]))
    text(plot_point_data_after(i),data_demod_nofalse(plot_point_data_after(i)),'o');
    title("�滻��Ĵ���λ�õĽ��������");
    legend('δ����������','û�д������','��������������');
    
    %% signal
    figure
    plot(signal_notcorrect,'b')
    hold on
    plot(signal_aftercorrect,'r+-')
    plot(signal_forcorrect,'g')
    axis(([plot_point_after(i)-50 plot_point_after(i)+50 -1.2 1.2]))
    text(plot_point_after(i),signal_notcorrect(plot_point_after(i)),'o');
    title("�滻��Ĵ���λ��");
    legend('δ����������','�����������','��������������');    
    %%
    figure
    plot(signal_notcorrect,'b')
    hold on
    plot(signal_forcorrect,'g')
    axis(([plot_point_after(i)-50 plot_point_after(i)+50 -1.2 1.2]))
    text(plot_point_after(i),signal_notcorrect(plot_point_after(i)),'o');
    title("�滻��Ĵ���λ��");
    legend('δ����������','��������������');  
    %%
    figure
    plot(signal_aftercorrect,'b')
    hold on
    plot(signal_forcorrect,'g')
    axis(([plot_point_after(i)-50 plot_point_after(i)+50 -1.2 1.2]))
    text(plot_point_after(i),signal_aftercorrect(plot_point_after(i)),'o');
    title("�滻��Ĵ���λ��");
    legend('�����������','��������������');  
    %%    
    figure
    plot(signal_notcorrect,'b')
    hold on
    plot(signal_aftercorrect,'r')
    axis(([plot_point_after(i)-50 plot_point_after(i)+50 -1.2 1.2]))
    text(plot_point_after(i),signal_notcorrect(plot_point_after(i)),'o');
    title("�滻��Ĵ���λ��");
    legend('δ����������','�����������');   
    %%
    figure
    plot(signal_nofalse,'b')
    hold on
    plot(signal_notcorrect,'r')
    plot(signal_aftercorrect,'g')
    axis(([plot_point_after(i)-50 plot_point_after(i)+50 -1.2 1.2]))
    text(plot_point_after(i),signal_notcorrect(plot_point_after(i)),'o');
    title("�滻��Ĵ���λ��");
    legend('δ�����ŵ�������','δ����������','�����������');    
    %%    
    figure
    plot(signal_nofalse,'b')
    hold on
    plot(signal_forcorrect,'r')
    axis(([plot_point_after(i)-50 plot_point_after(i)+50 -1.2 1.2]))
    text(plot_point_after(i),signal_forcorrect(plot_point_after(i)),'o');
    title("�滻��Ĵ���λ��");
    legend('δ�����ŵ�������','��������������');    
    %% 
    figure
    plot(signal_nofalse,'b')
    hold on
    plot(signal_notcorrect,'r')
    plot(signal_inf,'g')
    axis(([plot_point_after(i)-50 plot_point_after(i)+50 -1.2 1.2]))
    text(plot_point_after(i),signal_notcorrect(plot_point_after(i)),'o');
    title("�滻��Ĵ���λ��");
    legend('δ�����ŵ�������','δ����������','�ο��ŵ�����');
    
end


% %% �����滻�źŵĴ���λ�õ�����
% for i = 1:length(plot_point_forcorrect)
%     close all    
%     %% data
%     figure
%     plot(data_demod_beforecorrect,'b')
%     hold on
%     plot(data_demod_aftercorrect,'r+-')
%     plot(data_demod_forcorrect,'g')
%     axis(([plot_point_data_forcorrect(i)-50 plot_point_data_forcorrect(i)+50 -0.5 3.5]))
%     text(plot_point_data_forcorrect(i),signal_notcorrect(plot_point_data_forcorrect(i)),'o');
%     title("�滻��Ĵ���λ�õĽ��������");
%     legend('δ����������','�����������','��������������');
%     %%
%     figure
%     plot(data_demod_beforecorrect,'r')
%     hold on
%     plot(data_demod_inf,'g')
%     plot(data_demod_nofalse,'b+-')
%     axis(([plot_point_data_forcorrect(i)-50 plot_point_data_forcorrect(i)+50 -0.5 3.5]))
%     text(plot_point_data_forcorrect(i),signal_notcorrect(plot_point_data_forcorrect(i)),'o');
%     title("�滻��Ĵ���λ�õĽ��������");
%     legend('δ����������','�ο��ŵ�����','û�д������');
%     %%
%     figure
%     plot(data_demod_beforecorrect,'b')
%     hold on
%     plot(data_demod_nofalse,'r+-')
%     plot(data_demod_forcorrect,'g')
%     axis(([plot_point_data_forcorrect(i)-50 plot_point_data_forcorrect(i)+50 -0.5 3.5]))
%     text(plot_point_data_forcorrect(i),signal_notcorrect(plot_point_data_forcorrect(i)),'o');
%     title("�滻��Ĵ���λ�õĽ��������");
%     legend('δ����������','û�д������','��������������');
%     
%     %% signal
%     figure
%     plot(signal_notcorrect,'b')
%     hold on
%     plot(signal_aftercorrect,'r+-')
%     plot(signal_forcorrect,'g')
%     axis(([plot_point_forcorrect(i)-50 plot_point_forcorrect(i)+50 -1.2 1.2]))
%     text(plot_point_forcorrect(i),signal_notcorrect(plot_point_forcorrect(i)),'o');
%     title("�滻��Ĵ���λ��");
%     legend('δ����������','�����������','��������������');    
%     %%    
%     figure
%     plot(signal_notcorrect,'b')
%     hold on
%     plot(signal_aftercorrect,'r')
%     axis(([plot_point_forcorrect(i)-50 plot_point_forcorrect(i)+50 -1.2 1.2]))
%     text(plot_point_forcorrect(i),signal_notcorrect(plot_point_forcorrect(i)),'o');
%     title("�滻��Ĵ���λ��");
%     legend('δ����������','�����������');   
%     %%
%     figure
%     plot(signal_nofalse,'b')
%     hold on
%     plot(signal_notcorrect,'r')
%     plot(signal_aftercorrect,'g')
%     axis(([plot_point_forcorrect(i)-50 plot_point_forcorrect(i)+50 -1.2 1.2]))
%     text(plot_point_forcorrect(i),signal_notcorrect(plot_point_forcorrect(i)),'o');
%     title("�滻��Ĵ���λ��");
%     legend('δ�����ŵ�������','δ����������','�����������');    
%     %%    
%     figure
%     plot(signal_nofalse,'b')
%     hold on
%     plot(signal_forcorrect,'r')
%     axis(([plot_point_forcorrect(i)-50 plot_point_forcorrect(i)+50 -1.2 1.2]))
%     text(plot_point_forcorrect(i),signal_forcorrect(plot_point_forcorrect(i)),'o');
%     title("�滻��Ĵ���λ��");
%     legend('δ�����ŵ�������','��������������');    
%     %% 
%     figure
%     plot(signal_nofalse,'b')
%     hold on
%     plot(signal_notcorrect,'r')
%     plot(signal_inf,'g')
%     axis(([plot_point_forcorrect(i)-50 plot_point_forcorrect(i)+50 -1.2 1.2]))
%     text(plot_point_forcorrect(i),signal_notcorrect(plot_point_forcorrect(i)),'o');
%     title("�滻��Ĵ���λ��");
%     legend('δ�����ŵ�������','δ����������','�ο��ŵ�����');
%     
% end
