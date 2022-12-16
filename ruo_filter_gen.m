function[new_origin_rate , new_d_rate , upf_transmit , dof_transmit , filter_transmit , upf_receive , dof_receive , filter_receive] ...
    = ruo_filter_gen(origin_rate,f_rate,times,filter_order,rp,rst)

    err = 1e20;
    for n = 1:100
        new_orignal_rate_tmp = (f_rate*n)/round(f_rate*n/origin_rate);
        llcm_transmit_tmp = f_rate*n;
        upf_transmit_tmp = llcm_transmit_tmp/new_orignal_rate_tmp;  

        new_d_rate_tmp = new_orignal_rate_tmp*times;
        for i = 1:times
            receive_rate = llcm_transmit_tmp*i;
            if mod(receive_rate,new_d_rate_tmp) == 0
                llcm_receive_tmp = receive_rate;
            end
        end
        upf_receive_tmp = llcm_receive_tmp/f_rate;

        err_tmp = abs(new_orignal_rate_tmp - origin_rate);
        if (err_tmp < err) && (upf_transmit_tmp <= 200) && (upf_receive_tmp <= 200)
            err = err_tmp; 
            new_origin_rate = new_orignal_rate_tmp;
            new_d_rate = new_d_rate_tmp;
            llcm_transmit = llcm_transmit_tmp;
            llcm_receive = llcm_receive_tmp;

            upf_transmit = upf_transmit_tmp;
            dof_transmit = llcm_transmit/f_rate;

            upf_receive = upf_receive_tmp;
            dof_receive = llcm_receive/new_d_rate;
            dof_receive = round(dof_receive);
        end
    end

    upf_transmit = round(upf_transmit);
    dof_transmit = round(dof_transmit);
    upf_receive = round(upf_receive);
    dof_receive = round(dof_receive);
    bw = new_origin_rate/2;
    ups_rate_transmit = new_origin_rate*upf_transmit;
    ups_rate_receive = f_rate*upf_receive;

    filter_transmit = firceqrip(filter_order,bw/(ups_rate_transmit/2),[rp rst],'passedge');
    filter_transmit = filter_transmit./norm(filter_transmit,2)*sqrt(bw/ups_rate_transmit*2); % filter used in function sam_rate_con in transmit part

    filter_receive = firceqrip(filter_order,bw/(ups_rate_receive/2),[rp rst],'passedge'); % Impulse function
    filter_receive = filter_receive./norm(filter_receive,2)*sqrt(bw/ups_rate_receive*2); % filter used in function sam_rate_con in receive part

end