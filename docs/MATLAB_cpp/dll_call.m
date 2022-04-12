loadlibrary('out\x64_Release\Wrapper\Wrapper.dll', 'include\BS_Simulink.h');

BS = calllib('Wrapper','new_BS');

calllib('Wrapper','BS_initialize', BS, 10);
expect45 = calllib('Wrapper','BS_sum', BS);

expect45

