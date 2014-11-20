% Compare results with known references
% Most of the time, it consists in running the content of the HOWTO files
% and compare the matrices with saved version previously computed.

clear all;


%  % HOWTO_egg --------------------------------------------------------------------
%  disp('Checking peakdet ...');
%  [egg, fs] = wavread('glottalized-m1.egg');
%  [results_matrix, SdSIG, SIG] = peakdet(egg, fs, 200, 3);
%  REF = load('test/references/peakdet.mat');
%  test_checkwithreference(results_matrix, REF.results_matrix);
%  
%  disp('Checking decom ...');
%  [egg, fs] = wavread(['howtos' filesep 'DB-cresc-a-D4-m1.egg']);
%  degg=diff(egg);
%  degg(length(egg))=degg(length(egg)-1);
%  degg = degg./max(abs(degg));
%  [decom_Oqeggmess, decom_f0mess, decom_npicferms, decom_npicouvers, decom_atimes] = decom_sig(degg, fs);
%  REF = load('test/references/decom.mat');
%  test_checkwithreference(decom_Oqeggmess, REF.decom_Oqeggmess);
%  test_checkwithreference(decom_f0mess, REF.decom_f0mess);
%  test_checkwithreference(decom_npicferms, REF.decom_npicferms);
%  test_checkwithreference(decom_npicouvers, REF.decom_npicouvers);
%  test_checkwithreference(decom_atimes, REF.decom_atimes);
%  disp('Checking decom howard ...');
%  [howard_Oqeggmess, howard_f0mess, howard_atimes] = oq_egg_sig(egg, fs, [], 'howard');
%  test_checkwithreference(howard_Oqeggmess, REF.howard_Oqeggmess);
%  test_checkwithreference(howard_f0mess, REF.howard_f0mess);
%  test_checkwithreference(howard_atimes, REF.howard_atimes);
%  
%  
%  % HOWTO_sinusoidal -------------------------------------------------------------

fname = '0011.arctic_bdl1';
[wav, fs] = wavread([fname '.wav']);
f0s = load([fname '.f0.txt']);

%  disp('Checking sin_analysis peak picking ...');
%  opt = sin_analysis();
%  opt.fharmonic  = false;
%  opt.use_ls     = false;
%  opt.resyn      = true; % Use the internal OLA method for the resynthesis
%  opt.debug      = 0;
%  [frames syn_sm] = sin_analysis(wav, fs, f0s, opt);
%  REF = load('test/references/sin_analysis-pp.mat');
%  test_checkwithreference(frames, REF.frames);
%  test_checkwithreference(syn_sm, REF.syn_sm);
%  
%  disp('Checking sin_analysis HM and Least Squares solution...');
%  opt = sin_analysis();
%  opt.fharmonic  = true;
%  opt.fadapted   = false;
%  opt.use_ls     = true;
%  opt.debug      = 0;
%  frames = sin_analysis(wav, fs, f0s, opt);
%  REF = load('test/references/sin_analysis-hm-ls.mat');
%  test_checkwithreference(frames, REF.frames);
%  hmsynopt = hm_synthesis4();
%  hmsynopt.debug = 0;
%  syn_hm = hm_synthesis4(frames, length(wav), fs, hmsynopt); % Use the harmonic resynthesis
%  test_checkwithreference(syn_hm, REF.syn_hm);
%  
%  disp('Checking ahm_air_analysis2 ...');
%  optair = ahm_air_analysis2();
%  optair.do_final_ahm_step = false;
%  optair.debug = 0;
%  [f0sair, frames] = ahm_air_analysis2(wav, fs, f0s, optair);
%  REF = load('test/references/ahm_air_analysis2.mat');
%  test_checkwithreference(f0sair, REF.f0sair);
%  test_checkwithreference(frames, REF.frames);


% HOWTO_HMPD -------------------------------------------------------------------
disp('Checking hmpd_analysis_harmonic ...');
hmpdopt = hmpd_analysis();
hmpdopt.f0min   = 75; % To adapt to the analyzed voice
hmpdopt.f0max   = 220;% To adapt to the analyzed voice
hmpdopt.debug = 0;
hmpdframes = hmpd_analysis_harmonic(wav, fs, f0s, hmpdopt);
REF = load('test/references/hmpd_analysis_harmonic.mat');
test_checkwithreference(hmpdframes, REF.hmpdframes);

disp('Checking hmpd_analysis_features uncompressed features ...');
[hmpdf0s, hmpdAE, hmpdPDM, hmpdPDD] = hmpd_analysis_features(hmpdframes, fs, hmpdopt);
REF = load('test/references/hmpd_analysis_features.mat');
test_checkwithreference(hmpdf0s, REF.hmpdf0s);
test_checkwithreference(hmpdAE, REF.hmpdAE);
test_checkwithreference(hmpdPDM, REF.hmpdPDM, false);
test_checkwithreference(hmpdPDD, REF.hmpdPDD, false);

disp('Checking hmpd_synthesis ...');
hmpdPDD = zeros(size(hmpdPDD)); % PDD synthesis cannot be tested numerically because it generates noise which will be different from one run to another.
synopt = hmpd_synthesis();
synopt.enc = hmpdopt; 
synopt.usemex = hmpdopt.usemex; % Speed up with mex function
synopt.debug = 0;
syn_pdd = hmpd_synthesis(hmpdf0s, hmpdAE, [], hmpdPDD, fs, length(wav), synopt);
REF = load('test/references/hmpd_synthesis.mat');
test_checkwithreference(syn_pdd, REF.syn_pdd, false);
syn_pdmpdd = hmpd_synthesis(hmpdf0s, hmpdAE, hmpdPDM, hmpdPDD, fs, length(wav), synopt);
test_checkwithreference(syn_pdmpdd, REF.syn_pdmpdd, false);

disp('Checking hmpd_analysis_features compressed features ...');
hmpdopt.amp_enc_method=2; hmpdopt.amp_log=true; hmpdopt.amp_order=39;
hmpdopt.pdd_log=true; hmpdopt.pdd_order=12;% MFCC-like phase variance
hmpdopt.pdm_log=true; hmpdopt.pdm_order=24;% Number of log-Harmonic coefs
[hmpdf0s, hmpdAE, hmpdPDM, hmpdPDD] = hmpd_analysis_features(hmpdframes, fs, hmpdopt);
REF = load('test/references/hmpd_analysis_features_compressed.mat');
test_checkwithreference(hmpdf0s, REF.hmpdf0s);
test_checkwithreference(hmpdAE, REF.hmpdAE);
test_checkwithreference(hmpdPDM, REF.hmpdPDM, false);
test_checkwithreference(hmpdPDD, REF.hmpdPDD, false);


% HOWTO_envelope ---------------------------------------------------------------
% TODO 

% HOWTO_formant ----------------------------------------------------------------
% TODO 

% HOWTO_glottalsource ----------------------------------------------------------
% TODO 

% HOWTO_spectra ----------------------------------------------------------------
% TODO 

disp(' ');
disp('Regression tests completed.');
disp('(If no warning appeared and mean relative errors are small, it means that all went well!)');
