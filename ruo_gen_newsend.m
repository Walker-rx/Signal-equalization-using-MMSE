function signal_new = ruo_gen_newsend(signal_origin,M) 
    signal_new = signal_origin;   
    if mod(M,2) == 0
        times = M/2 ;
        for i = 1:times-1
            signal_new( 2*(i-1) <= signal_new & signal_new < 2*i ) = 2*i-1;
            signal_new( -2*i <= signal_new & signal_new < -2*(i-1) ) = -2*i+1;
        end
        signal_new( 2*(times-1) <= signal_new  ) = 2*times-1;
        signal_new(  signal_new < -2*(times-1) ) = -2*times+1;
    else
        times = floor(M/2);
        for j = 0:times-1
            signal_new( 2*j-1 <= signal_new & signal_new < 2*j+1 ) = 2*j;
            signal_new( -(2*j+1) <= signal_new & signal_new < -(2*j-1) ) = -2*j;
        end
        signal_new( 2*times-1 <= signal_new  ) = 2*times;
        signal_new(  signal_new < -(2*times-1) ) = -2*times;
    end
end