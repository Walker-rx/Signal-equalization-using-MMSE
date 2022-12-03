function [length_loop,ps_loop,pn_loop,errornum_ls_loop,error_location,data_demod_ls] = ruo_calculate_ser(data,signal_ori,signal_received,pilot_length,zero_length,data_length ...
                                                 ,fin_syn_point,times,ls_order)
                                             
    signal_downsample = signal_received(fin_syn_point:times:end);
    noise = signal_downsample(pilot_length+100:pilot_length+100+zero_length/2-1);
    pn_loop = bandpower(noise);
    data_received = signal_downsample(pilot_length+zero_length+1:pilot_length+zero_length+1+data_length-1);
    p = bandpower(data_received);
    ps_loop = p - pn_loop;
    
    signal_demod_ls = ruo_signal_equal_ls(signal_ori,signal_received,times,fin_syn_point,pilot_length,zero_length,ls_order);
    data_demod_ls = signal_demod_ls(pilot_length+zero_length+1:end);
    
    if length(data_demod_ls) < data_length
        compare_length = length(data_demod_ls);
        errornum_ls_loop = sum(data_demod_ls ~= data(1:compare_length));
        error_location = find(data_demod_ls ~= data(1:compare_length));
    else
        compare_length = data_length;
        errornum_ls_loop = sum(data_demod_ls(1:compare_length) ~= data);
        error_location = find(data_demod_ls(1:compare_length) ~= data);
    end
    length_loop = compare_length;
                
end