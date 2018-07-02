  %statistics of short cycles
   Ng_mean_VN = mean(Ng_per_p_);
   Ng2_mean_VN = mean(Ng2_per_p_);
   Ng4_mean_VN = mean(Ng4_per_p_);
   Ng_mean_CN = mean(Ng_per_u_);
   Ng2_mean_CN = mean(Ng2_per_u_);
   Ng4_mean_CN = mean(Ng4_per_u_);
   Ng_std_VN_= std(Ng_per_p_);
   Ng2_std_VN = std(Ng2_per_p_);
   Ng4_std_VN = std(Ng4_per_p_);
   Ng_std_CN=std(Ng_per_u_);
   Ng2_std_CN=std(Ng2_per_u_);
   Ng4_std_CN=std(Ng4_per_u_);

   %QC_LDPC_statistics
   Circulant_size=61;
   columns=10;
   rows=3;
   for i = 1:columns
   Ng_mean_VN_QC (i)= Ng_per_p_((i-1)*Circulant_size+1);
   end

  for i = 1:rows
   Ng_mean_CN_QC (i)= Ng_per_u_((i-1)*Circulant_size+1);
   end

  for i = 1:columns
   Ng2_mean_VN_QC (i)= Ng2_per_p_((i-1)*Circulant_size+1);
   end

  for i = 1:rows
   Ng2_mean_CN_QC (i)= Ng2_per_u_((i-1)*Circulant_size+1);
   end
  for i = 1:columns
   Ng4_mean_VN_QC (i)= Ng4_per_p_((i-1)*Circulant_size+1);
   end

  for i = 1:rows
   Ng4_mean_CN_QC (i)= Ng4_per_u_((i-1)*Circulant_size+1);
   end