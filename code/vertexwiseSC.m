function vertexwiseSC(study, subid, mm, n, stage)
% VERTEX-WISE STRUCTURAL COVARIANCE
%   This script assumes that the datasets are on your path.
%
%   Options:
%    - study: str ; which study do we want to analyze (i.e. {Insight, NUSDAST})
%    - subid: int ; row number (i.e. from 1:length of dataset)
%    - mm: int ; spatial smoothing in mm (i.e. one of {10, 20, 40})
%    - n: int ; index of which network to process (i.e. from 1:7)
%    - stage: int ; dictates what processing stage we're in (i.e. from 1:3)

  %% House-keeping
  switch stage
    case 1
      SC_preproc(study, mm)
    case 2
      SC_map(study, subid, mm, n)
    case 3
      SC_reduce(study, mm, n)
  end % endswitch
end % endfunction

function SC_preproc(study, mm)
  %% LOAD DATA (35 seconds)
  switch study
      case 'Insight'
          load InsightCT_0mm.mat
          load InsightBehData.mat
      case 'NUSDAST'
          load NUSDAST_CT_0mm.mat
          load NUSDASTBehData.mat
  end%switch
  load CIVETavsurf.mat
  load CIVETmask.mat
  load Yeo7networks_info.mat

  %% ANALYSIS
  residvars = zeros(size(Y, 1),3);
  residvars(:, 1) = beh.Age;
  sex = double(term(beh.Sex));
  sex = sex(:, 1);
  sex(sex==0)=2;
  residvars(:, 2) = sex;
  r = struct();
  
  groups = unique(beh.Group);
  Group = term(beh.Group);
  ngroups = size(groups,1);
  
  % Smoothing - 1.2 seconds (1mm) to 230 seconds (40mm)
  smoothed = SurfStatSmooth(Y, avsurf, mm/2);
  for n = 1:size(info.abbreviation,1)
      display(info.abbreviation(n))
      nverts = info.numverts(n);
      
      % Parcellation (< 1 sec)
      parcellated_network = smoothed(:, info.ROIverts==n);
      parcellated_datamask = mask(:, info.ROIverts==n);
      network_meanthick = mean(double(parcellated_network(:, parcellated_datamask)),2);
      residvars(:,3) = network_meanthick;
  
      % Covariates (< 1 sec)
      M = 1 + term(residvars);
      slm = SurfStatLinMod(parcellated_network, M);
      resid = parcellated_network - slm.X*slm.coef;
  
      outfname = "./data/resid/vertexConnectivity_" + study + "_" + mm + "_" + info.abbreviation(n) + "_resid.mat";
      save(outfname, 'resid', '-v7.3')
  
  end%for

  % Saved above: residuals for the network
  % Previously saved as:
  %r.(sname).Y = smoothed;
  %r.(sname).(netname).Yn = parcellated_network;
  %r.(sname).(netname).Yn_meanCT = meanthick;
  %r.(sname).(netname).Yn_resid = resid;

end%function

function SC_map(study, subid, mm, n)
% Computes features for a subject

  switch study
      case 'Insight'
          load InsightBehData.mat
      case 'NUSDAST'
          load NUSDASTBehData.mat
  end%switch
  load Yeo7networks_info.mat
  nverts = info.numverts(n);
  
  Group = term(beh.Group);
  groups = unique(beh.Group);
  ngroups = size(groups,1);

  infname = "./data/resid/vertexConnectivity_" + study + "_" + mm + "_" + info.abbreviation(n) + "_resid.mat";
  load(infname)

  % Group Matrices (~ 25 seconds)
  corr = zeros(nverts, nverts, ngroups);
  for gr = 1:ngroups 
      grname = groups(gr);
      grCT = resid(Group.(grname) == 1, :);
      % Correlations (~ 9 - 20 seconds)
      corr(:,:,gr) = corrcoef(grCT);
      corr(:,:,gr) = corr(:,:,gr) - diag(diag(corr(:,:,gr)));
  end%for

  resid_LOO = resid;
  resid_LOO(subid, :) = [];
  corrLOO = corrcoef(resid_LOO);

  if strcmp(beh.Group{subid}, 'Control')
    W = corr(:,:,1) - corrLOO;
  elseif strcmp(beh.Group{subid}, 'Patient')
    W = corr(:,:,2) - corrLOO;
  end
  normW = W - min(W(:));
  normW = normW./max(normW(:));          

  gr = struct();
  % Modularity (~30 mins)
  display("Modularity")
  [gr.Ci, gr.Q] = modularity_und(normW);
  display("Module Degree")
  gr.Z = module_degree_zscore(normW, gr.Ci);
  display("Participation Coefficient")
  gr.PC = participation_coef(normW,gr.Ci);
  display("Strength")
  gr.S = strengths_und(normW);
  display("Rich Club Coefficients (Weighted)")
  gr.RCW = rich_club_wu(normW);

  outfname = "./data/feat/vertexConnectivity_" + study + "_" + mm + "_" + info.abbreviation(n) + "_" + subid + "_features.mat";
  save(outfname, 'gr', '-v7.3')

  % Saved above: graph network features
  % Previously saved as:
  %r.(sname).(netname).SC_control = corr(:,:,1);
  %r.(sname).(netname).SC_patient = corr(:,:,2);
  %r.(sname).(netname).graph = gr;

end%function

function SC_reduce(study, mm, n)
% Merges feature tables and performs group comparisons
  switch study
      case 'Insight'
          load InsightBehData.mat
      case 'NUSDAST'
          load NUSDASTBehData.mat
  end%switch
  load Yeo7networks_info.mat
  nverts = info.numverts(n);
  nsubs = size(beh.Group, 1);
  
  Group = term(beh.Group);
  groups = unique(beh.Group);
  ngroups = size(groups,1);

  g = struct();
  g.Ci = zeros(nverts, nsubs);
  g.Q = zeros(nsubs, 1);
  g.Z = zeros(nverts, nsubs);
  g.PC = zeros(nverts, nsubs);
  g.S = zeros(nverts, nsubs);
  g.RCW = zeros(nverts, nsubs);

  infname_base = "./data/feat/vertexConnectivity_" + study + "_" + mm + "_" + info.abbreviation(n) + "_";
  for subid = 1:nsubs
    infname = infname_base + subid + "_features.mat";
    load(infname)

    g.Ci(:, subid) = gr.Ci;
    g.Q(subid) = gr.Q;
    g.Z(:, subid) = gr.Z;
    g.PC(:, subid) = gr.PC;
    g.S(:, subid) = gr.S;
    g.RCW(:, subid) = gr.RCW;
  end%for
  gr = g;

  % Compare Groups
  t.mod.patient = gr.Q(logical(Group.Patient));
  t.mod.control = gr.Q(logical(Group.Control));
  
  [t.mod.h, t.mod.p, t.mod.ci, t.mod.stats] = ttest2(t.mod.control, t.mod.patient);
  [t.mod.hu, t.mod.pu] = ttest2(t.mod.control, t.mod.patient,'Vartype','unequal');

  % Correlate with Symptoms/Cognition
  % TODO: add all other graph features alongside Modularity
  c.mod.labels = {'Modularity', 'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention', 'Social Cognition'};
  [c.mod.r,c.mod.p] = corrcoef([gr.Q, beh.VerbMem, beh.VisMem, beh.WorkMem, beh.ProcSpeed, beh.ExecFunc, beh.Att, beh.SocCog], 'Rows', 'pairwise');

  c.mod.labels_clin = {'Modularity', 'SAPS', 'SANS'};
  pt = logical(Group.Patient);
  [c.mod.r_clin,c.mod.p_clin] = corrcoef([gr.Q(pt), beh.SAPS_35(pt), beh.SANS_27(pt)], 'Rows', 'pairwise');

  % Data organization
  sname = ['smooth' num2str(mm)];
  netname = info.abbreviation(n);
  netname = netname{1};
  r = struct();
  r.ttest = t;
  r.corr = c;

  outfname = "./data/stats/vertexConnectivity_" + study + "_" + mm + "_" + info.abbreviation(n) + "_results.mat";
  save(outfname, 'r', '-v7.3')

  % Saved above: t-test results and correlations
  % Previously saved as:
  %r.(sname).(netname).ttest = t;
  %r.(sname).(netname).corr = c;

end%function


