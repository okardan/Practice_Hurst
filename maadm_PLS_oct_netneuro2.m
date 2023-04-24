% Reference:
% Kardan, O., et al. (2023) Improvements in task performance after practice
% are associated with scale-free dynamics of brain activity.
% Network Neuroscience, 1-63. https://doi.org/10.1162/netn_a_00319

% Data used in this script are here: https://osf.io/zsxfj

% This script runs PLS analysis on CAST data parcel-wise H (Study 3)
% Results are used in Figure 5 of the paper and corresponding supplementary
% analyses using Craddock 392 parcellation, task-timings-regressed, or WLMF

clear all

addpath(genpath('~\NIFTI_tools')); % Copyright (c) 2009, Jimmy Shen see license.txt inside NIFTI_tools
shen268 = load_untouch_nii('~\NUBE_light\shen_1mm_268_parcellation.nii'); %Shen et al.,(2013) parcellation
cc392 = load_untouch_nii('~\NUBE_light\CC400.nii'); %Craddock et al.,(2012) parcellation

load('~\Shen_Parcel_H_data_studies1&3\MaadOct2023.mat');
incsubs =1:44;
rng('default');

%*****
cast1 =[]; cast2=[]; cast3=[]; cast4=[]; cast5=[]; cast6=[];
for i=1:44
    cast1 = [cast1; MaadOct2023{i, 1}.run1.Aprm];
    cast2 = [cast2; MaadOct2023{i, 1}.run2.Aprm];
    cast3 = [cast3; MaadOct2023{i, 1}.run3.Aprm];
    cast4 = [cast4; MaadOct2023{i, 1}.run4.Aprm];
    cast5 = [cast5; MaadOct2023{i, 1}.run5.Aprm];
    cast6 = [cast6; MaadOct2023{i, 1}.run6.Aprm];
end
Mc = mean([cast1 cast2 cast3 cast4 cast5 cast6],2);
% option.stacked_behavdata = [cast1 - Mc; cast2 - Mc; cast3 - Mc; cast4 - Mc; cast5 - Mc; cast6 - Mc];
avH = [];
for i =1:44
    avH = [avH;  nanmean([MaadOct2023{i, 1}.run1.shH(:,1) MaadOct2023{i, 1}.run2.shH(:,1) MaadOct2023{i, 1}.run3.shH(:,1),...
        MaadOct2023{i, 1}.run4.shH(:,1) MaadOct2023{i, 1}.run5.shH(:,1) MaadOct2023{i, 1}.run6.shH(:,1)])];
end
Mh = mean(avH,2);
%*******

addpath(genpath('~\PlS')); % download PLS containing plscmd and plsgui folders from https://github.com/McIntosh-Lab/PLS
MaadDiff =[]; MADperfBase=[];
for i=1:44
    MADperfBase = [MADperfBase; MaadOct2023{i, 1}.run1.Aprm];
    MaadDiff = [MaadDiff;MaadOct2023{i, 1}.run6.Aprm - MaadOct2023{i, 1}.run1.Aprm];
end
yy = MaadDiff; XX = [ones(44,1) MADperfBase];
[b,bint,classlab,rint,stats]  = regress(yy,XX);

groups{1} = [1:44];
nRuns = 6;Runs =[1:6];

% switch between Shen parcels or Craddock 392 (CC400) parcels
nParcels = 268; badParcels =[4,131,137,189,239,252]; goodps = setdiff(1:268,badParcels); shencradd = 'Shen'; %parcels with NaN values in Study 1 and 3 are excluded
% nParcels = 392; badParcels=[48]; goodps = setdiff(1:392,badParcels); shencradd = 'Cradd';


%%
datamat_lst = cell(length(groups),1);
%###########********************
delaprimes = classlab; % or raw yy %adjusted change in perf or non-adjusted for supplementary
delaprimes1=[];delaprimes2=[];

for g = 1:numel(groups)
    for c = 1:nRuns
        for s = 1:numel(groups{g})
            ss = incsubs(s);
            
            %             hmat1 = MaadOct2023{ss, 1}.run1.craddH(goodps,1)'; %for supplementary Craddock analysis
            %             hmat2 = MaadOct2023{ss, 1}.run2.craddH(goodps,1)';
            %             hmat3 = MaadOct2023{ss, 1}.run3.craddH(goodps,1)';
            %             hmat4 = MaadOct2023{ss, 1}.run4.craddH(goodps,1)';
            %             hmat5 = MaadOct2023{ss, 1}.run5.craddH(goodps,1)';
            %             hmat6 = MaadOct2023{ss, 1}.run6.craddH(goodps,1)';
            
            %             hmat1 = MaadOct2023{ss, 1}.run1.tskreggedshenH(goodps,1)'; %for supplementary task_regressed analysis
            %             hmat2 = MaadOct2023{ss, 1}.run2.tskreggedshenH(goodps,1)';
            %             hmat3 = MaadOct2023{ss, 1}.run3.tskreggedshenH(goodps,1)';
            %             hmat4 = MaadOct2023{ss, 1}.run4.tskreggedshenH(goodps,1)';
            %             hmat5 = MaadOct2023{ss, 1}.run5.tskreggedshenH(goodps,1)';
            %             hmat6 = MaadOct2023{ss, 1}.run6.tskreggedshenH(goodps,1)';
            
            
            %             hmat1 = MaadOct2023{ss, 1}.run1.shH(goodps,2)'; %c1 for supplementary WLMF analysis
            %             hmat2 = MaadOct2023{ss, 1}.run2.shH(goodps,2)';
            %             hmat3 = MaadOct2023{ss, 1}.run3.shH(goodps,2)';
            %             hmat4 = MaadOct2023{ss, 1}.run4.shH(goodps,2)';
            %             hmat5 = MaadOct2023{ss, 1}.run5.shH(goodps,2)';
            %             hmat6 = MaadOct2023{ss, 1}.run6.shH(goodps,2)';
            
            hmat1 = MaadOct2023{ss, 1}.run1.shH(goodps,1)';  % main analyses
            hmat2 = MaadOct2023{ss, 1}.run2.shH(goodps,1)';
            hmat3 = MaadOct2023{ss, 1}.run3.shH(goodps,1)';
            hmat4 = MaadOct2023{ss, 1}.run4.shH(goodps,1)';
            hmat5 = MaadOct2023{ss, 1}.run5.shH(goodps,1)';
            hmat6 = MaadOct2023{ss, 1}.run6.shH(goodps,1)';
            
            if c ==1
                vec = hmat1;
                %                 vec = hmat1 - Mh(s); mean-centered H for Revision#2
            end
            if c ==2
                vec = hmat2;
                %                 vec = hmat2- Mh(s);
            end
            if c ==3
                vec = hmat3;
                %                 vec = hmat3- Mh(s);
            end
            if c ==4
                vec = hmat4;
                %                 vec = hmat4- Mh(s);
            end
            if c ==5
                vec = hmat5;
                %                 vec = hmat5- Mh(s);
            end
            if c ==6
                vec = hmat6;
                %                 vec = hmat6- Mh(s);
            end
            datamat_lst{g} = [datamat_lst{g}; vec];
            s
        end
    end
end


num_subj = [length(groups{1})];
num_cond = nRuns;

option.method = 3; %Task PLS1  or behavioral PLS3
option.num_boot = 5000;
option.num_perm = 1000;
option.meancentering_type=[2]; % 2 is for only grand mean removed
option.cormode = 0; %	0. Pearson correlation 2. covaraince 4. cosine angle 6. dot product
option.stacked_behavdata = repmat(delaprimes,6,1);

result = pls_analysis(datamat_lst, num_subj, num_cond, option);

%put in the datamat_lst into results.  No longer automatic in new PLS
%version

result.datamat_lst = datamat_lst;
ps = result.perm_result.sprob
brainPs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)>3)
goodps(find(result.boot_result.compare_u(:,1)>3))
brainNs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)<-3)
goodps(find(result.boot_result.compare_u(:,1)<-3))
figure;
for k=1:nRuns
    scatter(k,result.boot_result.ulcorr_adj(k,1));hold on
    scatter(k,result.boot_result.orig_corr(k,1));
    scatter(k,result.boot_result.llcorr_adj(k,1));
end
hold off
%%
%save('CASTadj_Shen_H_behpls_result.mat','result');

%% beh PLS figure making
%load('CASTadj_Shen_H_behpls_result.mat')
for lv=1:1
    %%
    figure;
    %flipped sign on both sides of LV for ease of interpretation
    line(ones(300,1),linspace(0,-result.boot_result.orig_corr(1,1),300),'Color','r','LineWidth',13);
    hold on
    for k=1:6
        
        col = 'k';
        xpos = k ;
        
        line(xpos*ones(300,1),linspace(0,-result.boot_result.orig_corr(k,lv),300),'Color',col,'LineWidth',13);
        scatter(xpos,-result.boot_result.orig_corr(k,lv),'+k');
        line(xpos*ones(300,1),linspace(-result.boot_result.ulcorr_adj(k,lv),-result.boot_result.llcorr_adj(k,lv),300));
    end
    
    sigfeats = find(abs(result.boot_result.compare_u(:,lv))>3); % flipped the sign above in result.boot so flip in brain fig too
    sigfeatsP = find(result.boot_result.compare_u(:,lv)>3);
    sigfeatsN = find(result.boot_result.compare_u(:,lv)<-3);
    z1percent = length(sigfeats==1);
    result.boot_result.compare_u(sigfeats,lv)
    cbcovs = result.s.^2./sum(result.s.^2);
    title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(ps(lv))]);
    ylim([-1,1.2]); xlim([.5,(nRuns +.5)]);
    ylabel('Correlation with adjusted (\Delta A'')');
    set(gca,'XTick',1:6,'Xticklabel',{'1st CAST run','2nd CAST run','3rd CAST run',...
        '4th CAST run','5th CAST run','6th CAST run'},...
        'XtickLabelRotation',45,'FontSize',14);
    hold off
    %%
    sigids = goodps(sigfeats);
    sigidsP = goodps(sigfeatsP);
    sigidsN = goodps(sigfeatsN);
    gen_parc_shen = zeros(61,72,56);
    gen_parc_gord = zeros(61,72,56);
    gen_parc_cc4 = zeros(61,72,56);
    for vx=3:61-3
        for vy=3:72-3
            for vz=3:56-3
                xmm = (vx-30.5)*3.25; vxshen = round(91-xmm);  vxcc4 = round(32-xmm/3);
                ymm = (vy-41.6)*3.25; vyshen = round(ymm+127);  vycc4 = round(ymm/3 +44);
                zmm = (vz-23.3)*3.5; vzshen = round(zmm+73);  vzcc4 = round(zmm/3 + 25);
                gen_parc_shen(vx,vy,vz) = shen268.img(vxshen,vyshen,vzshen);
                gen_parc_cc4(vx,vy,vz) = cc392.img(vxcc4,vycc4,vzcc4);
            end
        end
    end
    gen_parc_shenSigs = gen_parc_shen;
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsP))=-3;%flipped sign on both sides of LV for ease of interpretation
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsN))=+3;
    gen_parc_shenSigs(~ismember(gen_parc_shen,sigids) & gen_parc_shen~=0)=0;
    nii = make_nii(gen_parc_shenSigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
    nii.hdr.dime.datatype =64; % used to make Figure 5 in the paper
    save_nii(nii,['~\behPLS_figures\deltaCASTregBase_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);
    
    % for the craddock version
    %     gen_parc_cc4Sigs = gen_parc_cc4;
    %     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsP))=+3;%
    %     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsN))=-3;
    %     gen_parc_cc4Sigs(~ismember(gen_parc_cc4,sigids) & gen_parc_cc4~=0)=0;
    %
    %    nii = make_nii(gen_parc_cc4Sigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
    %     nii.hdr.dime.datatype =64;
    %     save_nii(nii,['~\behPLS_figures\deltaCASTregBase_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);
    
end