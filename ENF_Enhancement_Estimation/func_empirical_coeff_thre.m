function thre = func_empirical_coeff_thre(signal_length,repeat_time)
temp_coeff = zeros(1,repeat_time);
for i = 1:repeat_time
   x             = randn(1,signal_length);
   y             = randn(1,signal_length);
   temp          = corrcoef(x,y);
   temp_coeff(i) = temp(1,2);
end
thre             = max(temp_coeff);
end


