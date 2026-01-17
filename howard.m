
% Ftemp is [N_a,N_z]
Ftemp_vec = reshape(Ftemp,[N_a*N_z,1]);




EVpre=reshape(VKron(aprimeindex,:),[N_aprimediff,N_a*N_z,N_z]); % last dimension is zprime
EVKrontemp=interp1(EVinterpindex1,EVpre,EVinterpindex2); % interpolate V as Policy points to the interpolated indexes
EVKrontemp=reshape(EVKrontemp,[N_aprime*N_a*N_z,N_z]);  % last dimension is zprime
EVKrontemp=EVKrontemp(tempmaxindex2,:);
EVKrontemp=EVKrontemp.*pi_z_howards;
EVKrontemp(isnan(EVKrontemp))=0;
EVKrontemp=reshape(sum(EVKrontemp,2),[N_a,N_z]);
VKron=Ftemp+DiscountFactorParamsVec*EVKrontemp;