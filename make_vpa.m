% Make and store VPA rectangular diffmat.

NN = (251:300);

for N = NN
    disp(N)
    [PD_vpa, D_vpa] = rectdiff_vpa(N);
    PD_vpa_store{N} = PD_vpa;
    D_vpa_store{N} = D_vpa;
end

save vpa_diffmats PD_vpa_store D_vpa_store