%VERTEX-WISE STRUCTURAL COVARIANCE
% Repeat once for each dataset (Insight, NUSDAST)

%% LOAD DATA (35 seconds)
study = 'Insight';
% study = 'NUSDAST';

switch study
    case 'Insight'
        load InsightCT_0mm.mat
        load InsightBehData.mat
    case 'NUSDAST'
        load NUSDAST_CT_0mm.mat
        load NUSDASTBehData.mat
end

load CIVETavsurf.mat
load CIVETmask.mat
load Yeo7networks_info.mat

%% ANALYSIS
scalespace = [1,5,10,15,25,30,35,40]; % in mm
residvars = zeros(size(Y,1),3);
residvars(:,1) = beh.Age;
sex = double(term(beh.Sex));
sex = sex(:,1);
sex(sex==0)=2;
residvars(:,2) = sex;
r = struct();

for mm = 1:size(scalespace,2) % maybe drop this?
    % Smoothing - 1.2 seconds (1mm) to 230 seconds (40mm)
    smoothed = SurfStatSmooth(Y, avsurf, scalespace(mm)/2);
    for n = 1:size(info.abbreviation,1)
        nverts = info.numverts(n);
        
        % Parcellation (< 1 sec)
        parcellated_network = smoothed(:,info.ROIverts==n);
        parcellated_datamask = mask(:, info.ROIverts==n);
        network_meanthick = mean(double(parcellated_network(:,parcellated_datamask)),2);
        residvars(:,3) = network_meanthick;

        % Covariates (< 1 sec)
        M = 1;
        for var = 1:size(residvars,2)
            M = M + term(residvars(var));
        end
        slm = SurfStatLinMod(parcellated_network,M);
        resid = parcellated_network - slm.X*slm.coef;

        % Group Matrices (~ 25 seconds)
        groups = unique(beh.Group);
        Group = term(beh.Group);
        ngroups = size(groups,1);
        corr = zeros(nverts, nverts, ngroups);
        for gr = 1:ngroups     
            grname = groups(gr);
            grCT = resid(Group.(grname) == 1, :);
            % Correlations (~ 9 - 20 seconds)
            corr(:,:,gr) = corrcoef(grCT);
            corr(:,:,gr) = corr(:,:,gr) - diag(diag(corr(:,:,gr)));
        end
        
        gr.Ci = zeros(nverts,size(grCT,1));
        gr.Q = zeros(size(grCT,1),1);
        gr.Z = zeros(nverts,size(grCT,1));
        gr.PC = zeros(nverts,size(grCT,1));
        gr.E = zeros(size(grCT,1),1);

        for subj = 1:size(grCT,1)    
            % Jackknife Leave One Out (~ 15 seconds)
            grCT_LOO = grCT;
            grCT_LOO(subj,:) = [];
            corrLOO = corrcoef(grCT_LOO);
            if strcmp(beh.Group{subj}, 'Control')
                W = corr(:,:,1) - corrLOO;
            elseif strcmp(beh.Group{subj}, 'Patient')
                W = corr(:,:,2) - corrLOO;
            end
            normW = W - min(W(:));
            normW = normW./max(normW(:));          
            
            % Modularity (~30 mins)
            [gr.Ci(:,subj),gr.Q(subj,1)] = modularity_und(normW);
            gr.Z(:,subj) = module_degree_zscore(normW,gr.Ci);
            gr.PC(:,subj) = participation_coef(normW,gr.Ci);
            
            % Global Efficiency (3 days)
            gr.E(subj,1) = efficiency(normW);
        end
        
        % Compare Groups
        t.mod.patient = gr.Q(logical(Group.Patient));
        t.mod.control = gr.Q(logical(Group.Control));
        t.eff.patient = gr.E(logical(Group.Patient));
        t.eff.ccontrol = gr.E(logical(Group.Control));
        [t.mod.h, t.mod.p, t.mod.ci, t.mod.stats] = ttest2(t.mod.control, t.mod.patient);
        [t.mod.hu, t.mod.pu] = ttest2(t.mod.control, t.mod.patient,'Vartype','unequal');
        [t.eff.h, t.eff.p, t.eff.ci, t.eff.stats] = ttest2(t.eff.control, t.eff.patient);
        [t.eff.hu, t.eff.pu] = ttest2(t.eff.control, t.eff.patient,'Vartype','unequal');
        
        % Correlate with Symptoms/Cognition
        c.labels = {'Modularity', 'Global Efficiency', 'SAPS', 'SANS', 'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention', 'Social Cognition'};
        [c.r,c.p] = corrcoef([gr.Q,gr.E],[beh.SAPS_35, beh.SANS_27, beh.VerbMem, beh.VisMem, beh.WorkMem, beh.ProcSpeed, beh.ExecFunc, beh.Att, beh.SocCog]);
        
        % Data organization
        sname = ['smooth' num2str(scalespace(mm))];
        netname = info.abbreviation(n);
        r.(sname).Y = smoothed;
        r.(sname).(netname).Yn = parcellated_network;
        r.(sname).(netname).Yn_meanCT = meanthick;
        r.(sname).(netname).Yn_resid = resid;
        r.(sname).(netname).SC_control = corr(:,:,1);
        r.(sname).(netname).SC_patient = corr(:,:,2);
        r.(sname).(netname).graph = gr;
        r.(sname).(netname).ttest = t;
        r.(sname).(netname).corr = c;
    end
end

save vertexConnectivity_results.mat r
