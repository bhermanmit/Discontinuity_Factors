function plot_all(phimat,phi_SA,phi_both,phi_both_DF,phi_SA_DFs)

% grap fluxes by group
phimat_SA = reshape(phi_SA,2,length(phi_SA)/2);
phimat_SA_DFs = reshape(phi_SA_DFs,2,length(phi_SA_DFs)/2);
phimat_both = reshape(phi_both,2,length(phi_SA)/2);
phimat_both_DF = reshape(phi_both_DF,2,length(phi_both_DF)/2);

% x
x_hom = linspace(-21.42,21.42,length(phi_SA)/2);
x_het = linspace(-21.42,21.42,68);

% renormalize fluxes
phimat_SA(1,:) = phimat_SA(1,:)/trapz(x_hom,phimat_SA(1,:));
phimat_SA(2,:) = phimat_SA(2,:)/trapz(x_hom,phimat_SA(2,:));
phimat_SA_DFs(1,:) = phimat_SA_DFs(1,:)/trapz(x_hom,phimat_SA_DFs(1,:));
phimat_SA_DFs(2,:) = phimat_SA_DFs(2,:)/trapz(x_hom,phimat_SA_DFs(2,:));
phimat_both(1,:) = phimat_both(1,:)/trapz(x_hom,phimat_both(1,:));
phimat_both(2,:) = phimat_both(2,:)/trapz(x_hom,phimat_both(2,:));
phimat_both_DF(1,:) = phimat_both_DF(1,:)/trapz(x_hom,phimat_both_DF(1,:));
phimat_both_DF(2,:) = phimat_both_DF(2,:)/trapz(x_hom,phimat_both_DF(2,:));
phimat{1}(1,:) = phimat{1}(1,:)/(sum(phimat{1}(1,:)) + sum(phimat{2}(1,:))); 
phimat{2}(1,:) = phimat{2}(1,:)/(sum(phimat{1}(1,:)) + sum(phimat{2}(1,:)));
phimat{1}(2,:) = phimat{1}(2,:)/(sum(phimat{1}(2,:)) + sum(phimat{2}(2,:)));
phimat{2}(2,:) = phimat{2}(2,:)/(sum(phimat{1}(2,:)) + sum(phimat{2}(2,:)));


% OpenMC heterogeneous flux (no need to divide by volume with uniform mesh)
phi_het_both = [
9.20350E-02
0.326649
9.24465E-02
0.326818
9.35383E-02
0.327014
9.49514E-02
0.327190
9.82106E-02
0.326870
9.96184E-02
0.327943
0.100455
0.329894
0.101120
0.332059
0.101332
0.335134
0.102824
0.337682
0.106958
0.338753
0.107940
0.342239
0.105931
0.348007
0.106152
0.353047
0.107570
0.357404
0.109706
0.361969
0.114837
0.364683
0.116473
0.370407
0.114884
0.378425
0.115494
0.385991
0.117499
0.392463
0.120085
0.399247
0.125924
0.404066
0.127860
0.412226
0.126041
0.422817
0.126638
0.432941
0.129072
0.441103
0.130288
0.451378
0.131144
0.461081
0.130406
0.473416
0.126519
0.486509
0.124166
0.501058
0.120885
0.514504
0.117472
0.531777
0.109654
0.551177
0.107067
0.567611
0.106569
0.580228
0.108230
0.592743
0.113048
0.601571
0.115353
0.612719
0.116067
0.623884
0.116987
0.635860
0.116811
0.647270
0.119486
0.656827
0.126772
0.662109
0.128046
0.671668
0.123210
0.684593
0.122771
0.695142
0.124084
0.702970
0.127247
0.709759
0.134999
0.711994
0.136202
0.719017
0.130849
0.730063
0.130093
0.738191
0.131148
0.743866
0.134167
0.748026
0.142132
0.747751
0.143024
0.752601
0.137247
0.761358
0.136427
0.766254
0.138698
0.767781
0.139440
0.770387
0.140099
0.772019
0.138380
0.775173
0.132635
0.780321
0.130120
0.783429
0.128487
0.784539
0.127715
0.785216
];

phi_het_SA16 = [
    0.297552
1.10838
0.298706
1.10755
0.301779
1.10647
0.305714
1.10372
0.314808
1.09812
0.318008
1.09601
0.318238
1.09669
0.318471
1.09703
0.316457
1.09886
0.318691
1.09692
0.328136
1.08947
0.328019
1.08957
0.318091
1.09686
0.315399
1.09927
0.315335
1.09941
0.317968
1.09731
0.327712
1.08979
0.327877
1.08980
0.318112
1.09733
0.315645
1.09969
0.315728
1.09985
0.318433
1.09773
0.328209
1.09044
0.328129
1.09044
0.318640
1.09790
0.316552
1.09981
0.318601
1.09759
0.318381
1.09711
0.318013
1.09641
0.314802
1.09860
0.305580
1.10432
0.301691
1.10739
0.298841
1.10855
0.297646
1.10933
];

phi_het_SA30 = [
    0.182014
1.09164
0.183019
1.09134
0.185476
1.09054
0.189438
1.08862
0.198110
1.08383
0.201178
1.08237
0.200732
1.08348
0.200550
1.08429
0.198067
1.08629
0.200339
1.08483
0.210016
1.07814
0.209631
1.07820
0.199330
1.08480
0.196345
1.08680
0.196237
1.08699
0.199012
1.08495
0.209086
1.07822
0.209051
1.07848
0.199026
1.08538
0.196213
1.08748
0.196249
1.08730
0.199283
1.08531
0.209486
1.07857
0.209832
1.07839
0.200293
1.08508
0.198042
1.08642
0.200667
1.08436
0.200979
1.08361
0.201234
1.08265
0.198021
1.08399
0.189336
1.08869
0.185518
1.09030
0.182965
1.09088
0.181842
1.09127
];

% Separate het fluxes
phimat_het_both = reshape(phi_het_both,2,68);
phimat_het_SA16 = reshape(phi_het_SA16,2,34);
phimat_het_SA30 = reshape(phi_het_SA30,2,34);

% Interpolate Het Flux to same grid
phimat_het_both_grid(1,:) = interp1(x_het,phimat_het_both(1,:),x_hom);
phimat_het_both_grid(2,:) = interp1(x_het,phimat_het_both(2,:),x_hom);
phimat_het_SA16_grid(1,:) = interp1(x_het(1:34),phimat_het_SA16(1,:),x_hom(1:length(x_hom)/2),'linear','extrap');
phimat_het_SA16_grid(2,:) = interp1(x_het(1:34),phimat_het_SA16(2,:),x_hom(1:length(x_hom)/2),'linear','extrap');
phimat_het_SA30_grid(1,:) = interp1(x_het(35:68),phimat_het_SA30(1,:),x_hom(length(x_hom)/2+1:length(x_hom)),'linear','extrap');
phimat_het_SA30_grid(2,:) = interp1(x_het(35:68),phimat_het_SA30(2,:),x_hom(length(x_hom)/2+1:length(x_hom)),'linear','extrap');

% Normalize het flux
phimat_het_both_grid(1,:) = phimat_het_both_grid(1,:)/trapz(x_hom,phimat_het_both_grid(1,:));
phimat_het_both_grid(2,:) = phimat_het_both_grid(2,:)/trapz(x_hom,phimat_het_both_grid(2,:));
phimat_het_SA16_grid(1,:) = phimat_het_SA16_grid(1,:)/trapz(x_hom,[phimat_het_SA16_grid(1,:),phimat_het_SA30_grid(1,:)]);
phimat_het_SA30_grid(1,:) = phimat_het_SA30_grid(1,:)/trapz(x_hom,[phimat_het_SA16_grid(1,:),phimat_het_SA30_grid(1,:)]);
phimat_het_SA16_grid(2,:) = phimat_het_SA16_grid(2,:)/trapz(x_hom,[phimat_het_SA16_grid(2,:),phimat_het_SA30_grid(2,:)]);
phimat_het_SA30_grid(2,:) = phimat_het_SA30_grid(2,:)/trapz(x_hom,[phimat_het_SA16_grid(2,:),phimat_het_SA30_grid(2,:)]);

% Reference divided by SA
ref_div_SA = phimat_het_both_grid./([phimat_het_SA16_grid,phimat_het_SA30_grid]);
ref_div_SA(1,:) = ref_div_SA(1,:)/trapz(x_hom,ref_div_SA(1,:));
ref_div_SA(2,:) = ref_div_SA(2,:)/trapz(x_hom,ref_div_SA(2,:));

% Synthesized Homogeneous SA Form Function (Watch difference group indices)
syn(2,:) = phimat_both_DF(2,:).*[phimat_het_SA16_grid(1,:),phimat_het_SA30_grid(1,:)];
syn(1,:) = phimat_both_DF(1,:).*[phimat_het_SA16_grid(2,:),phimat_het_SA30_grid(2,:)];
syn(2,:) = syn(2,:)/trapz(x_hom,syn(2,:));
syn(1,:) = syn(1,:)/trapz(x_hom,syn(1,:));

% plot Group 1 Flux
figure(1)
hold on
plot(x_hom,phimat_SA(1,:),'b-','LineWidth',2);
plot(x_hom,phimat_SA_DFs(1,:),'c-','LineWidth',2);
plot(x_hom,phimat_both(1,:),'r-','LineWidth',2);
plot(x_hom,phimat_both_DF(1,:),'k-','LineWidth',2);
plot(x_hom,phimat_het_both_grid(2,:),'g--','LineWidth',2);
legend('Single Assembly Homogenization k = 1.1645',...
       'Single Assembly Homogenization w/ ADFs k = 1.1425',...
       'Two Assembly Homogenization k = 1.1640',...
       'Two Assembly Homogenization w/ RDFs k = 1.1393',...
       'OpenMC Heterogeneous Flux k = 1.1393');
xlabel('x - distance [cm]');
ylabel('Flux [-]');
title('Homogeneous Group 1 Flux Comparison');


% plot Group 2 Flux
figure(2)
hold on
plot(x_hom,phimat_SA(2,:),'b-','LineWidth',2);
plot(x_hom,phimat_SA_DFs(2,:),'c-','LineWidth',2);
plot(x_hom,phimat_both(2,:),'r-','LineWidth',2);
plot(x_hom,phimat_both_DF(2,:),'k-','LineWidth',2);
plot(x_hom,phimat_het_both_grid(1,:),'g--','LineWidth',2);
legend('Single Assembly Homogenization k = 1.1645',...
       'Single Assembly Homogenization w/ ADFs k = 1.1425',...
       'Two Assembly Homogenization k = 1.1640',...
       'Two Assembly Homogenization w/ RDFs k = 1.1393',...
       'OpenMC Heterogeneous Flux k = 1.1393');
xlabel('x - distance [cm]');
ylabel('Flux [-]');
title('Homogeneous Group 2 Flux Comparison');

% plot thermal flux profile from two single assemblies
figure(3)
hold on
plot(x_hom(1:length(x_hom)/2),phimat_het_SA16_grid(1,:),'b-','LineWidth',2);
plot(x_hom(length(x_hom)/2+1:length(x_hom)),phimat_het_SA30_grid(1,:),'r-','LineWidth',2);
legend('UOX 1.6%','UOX 3.37%');
xlabel('x - distance [cm]');
ylabel('Thermal Flux [-]');
title('Thermal Flux (Form Functions) from Single Assembly Calculations');

% plot OpenMC two-node solution divided by SA FF Group 1
figure(4)
hold on
plot(x_hom,ref_div_SA,'LineWidth',2);
legend('Group 1','Group 1');
xlabel('x - distance [cm]');
ylabel('Ratio of Flux to Form Function [-]');
title('OpenMC two-node solution divided by Single Assembly Form Function');

% plot Homogenized w/ RDFs two node homogenous solution multiplied by SA
% Form Function
figure(5)
hold on
plot(x_hom,syn(1,:),x_hom,syn(2,:),'LineWidth',2);
legend('Group 1','Group 2');
xlabel('x - distance [cm]');
ylabel('Flux mulitplied by Form Function [-]')
title('Synthesized Homogeneous w/RDFs and Single Assembly Form Functions')

% plot Heterogeneous reference solution
figure(6)
hold on
plot(x_hom,phimat_het_both_grid(2,:),x_hom,phimat_het_both_grid(1,:),'LineWidth',2);
legend('Group 2','Group 1');
xlabel('x - distance [cm]');
ylabel('Flux [-]');
title('OpenMC Reference Heterogeneous Flux');

% plot Synthesized with Reference Comparison
figure(7)
hold on
plot(x_hom,phimat_het_both_grid(2,:),'k-','LineWidth',2)
plot(x_hom,syn(1,:),'r-','LineWidth',2)
legend('OpenMC Reference','Synthesized');
xlabel('x - distance [cm]');
ylabel('Flux [-]');
title('Comparison of Group 1 OpenMC Reference to Homogeneous Synthesized with SA Form Function');

% plot Synthesized with Reference Comparison
figure(8)
hold on
plot(x_hom,phimat_het_both_grid(1,:),'k-','LineWidth',2)
plot(x_hom,syn(2,:),'r-','LineWidth',2)
legend('OpenMC Reference','Synthesized');
xlabel('x - distance [cm]');
ylabel('Flux [-]');
title('Comparison of Group 2 OpenMC Reference to Homogeneous Synthesized with SA Form Function');