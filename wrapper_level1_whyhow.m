function allinput = wrapper_level1_whyhow(varargin)
% WRAPPER_LEVEL1_WHYHOW
%
%   USAGE: allinput = wrapper_level1_whyhow(varargin)
% __________________________________________________________________________
%  OUTPUT
%   allinput: cell array containing matlabbatch structures for running model
%   specification, estimation, and then calculating contrasts
% __________________________________________________________________________
%  OPTIONAL INPUTS (VARARGIN)
%   These are entered as `'name', value` argument pairs. Matching is not
%   case-sensitive and partial name matches are OK. To see default values,
%   run `wrapper_level1_whyhow` without an output variable.
%
% THESE VARARGINS ARE USED FOR FINDING THE NECESSARY FILES AND DIRECTORIES
% | :------------- | ------------------------------------------------------------------------------------------ |
% | studydir       | full path to directory containing subject data folders                                     |
% | behavpat       | search pattern for finding behavioral data within subject dirs, e.g., 'behav/whyhow*mat'   |
% | epipat         | search pattern for finding functional data file(s) within each run, e.g., 'sw*nii'         |
% | nuisancepat    | search pattern for finding nuisance regressor file within each run, e.g., 'rp*txt'         |
% | runpat         | searc h pattern for finding run directories within subject dirs, e.g.,'raw/BOLD_WhyHow*'    |
% | subpat         | search pattern for finding subject directories within studydir, e.g., 'Subject*'           |
% | brainmask      | full path to brain mask to use (leave empty for none)                                      |
%
% THESE VARARGINS ARE USED TO SPECIFY RELEVANT DETAILS ABOUT THE IMAGE AND BEHAVIORAL DATA BEING MODELED
% | :------------- | ------------------------------------------------------------------------------------------  |
% | is4D           | flag for 4D image file (0=No, 1=Yes)                                                        |
% | nskip          | number of initial TRs that have been removed (for adjusting stimulus onsets)                |
% | TR             | acquisition repetition time (in seconds)                                                    |
% | yesnokeys      | keys corresponding to yes/no responses (e.g., [1 2])                                        |
%
% THESE VARARGINS ARE USED TO CONFIGURE THE MODEL AND ESTIMATION METHODS
% | :------------- | ------------------------------------------------------------------------------------------  |
% | basename       | base name for the analysis (e.g., 'WhyHow_SmoothedData')                                    |
% | model          | string specifying model to use ('2x2' for full design, '1x2' to collapse faces/hands factor |
% | incl_err       | flag to include parametric covariate modeling blockwise variation in # of errors            |
% | incl_rt        | flag to include parametric covariate modeling blockwise variation in response time          |
% | armethod       | autocorrelation removal method (0=None, 1=AR(1), 2=Weighted Least Squares (WLS), 3=FAST)    |
% | HPF            | high-pass filter cutoff to use (in seconds)                                                 |
% | maskthresh     | implicit masking threshold (proportion of globals), default = 0.8                           |
% | fcontrast      | flag to compute omnibus F-contrast (useful for feature selection, e.g., in PPI analysis)    |
% | run_it_now     | flag to run analysis now (0=No, 1=Yes)                                                      |
%

% | SET DEFAULTS AND PARSE VARARGIN
% | ===========================================================================
defaults = {
'armethod',         2,                  ...
'basename',         'WHYHOWLOC',        ...
'behavpat',         'behav/whyhow*mat', ...
'brainmask',        '',                 ...
'epipat',           'file*nii',        ...
'fcontrast',        1,                  ...
'HPF',              100,                ...
'incl_err',         1,                  ...
'incl_rt',          1,                  ...
'is4D',             1,                  ...
'maskthresh',       -Inf,                ...
'model',            '2x2',              ...
'nskip',            0,                  ...
'nuisancepat',      'rp*txt',           ...
'runpat',           'raw/bol*',        ...
'studydir',         '/Users/bobspunt/Documents/fmri/keren',                ...
'subpat',           'sub*',         ...
'TR',               1.5,                  ...
'run_it_now',       0,                  ...
'yesnokeys',        [1 2]               ...
 };
vals = setargs(defaults, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals);


% | PATHS
% | ===========================================================================
[subdir, subnam] = files([studydir filesep subpat]);

% | ANALYSIS NAME
% | ===========================================================================
covidx          = find([incl_rt incl_err]);
armethodlabels  = {'NoAR1' 'AR1' 'WLS'};
covnames        = {'Duration' 'Errors' 'FoilErrors'};
if ~isempty(covidx)
    pmnames         = regexprep(covnames(covidx), '_', '');
    pmstr           = sprintf(repmat('_%s', 1, length(pmnames)), pmnames{:}); pmstr(1)= [];
else
    pmstr = 'None';
end
analysisname  = sprintf('%s_%s_Pmodby_%s_%s_%ds_%s', basename, model, ...
                        pmstr, armethodlabels{armethod + 1}, HPF, bspm_timestamp);
printmsg(analysisname, 'msgtitle', 'Analysis Name');

% | IMAGING PARAMETERS
% | ========================================================================
adjons          = TR*nskip;

% | SUBJECT LOOP
% | ===========================================================================
allinput = [];
for s = 1:length(subdir)

    % | Check Subject and Define Folders
    % | ========================================================================
    rundir      = files([subdir{s} filesep runpat]);
    if isempty(rundir), printmsg('Valid run directory not found, moving on...', 'msgtitle', subnam{s}); continue; end
    analysisdir = fullfile(subdir{s}, 'analysis', analysisname);
    if any([exist(fullfile(analysisdir, 'mask.img'), 'file') exist(fullfile(analysisdir, 'mask.nii'), 'file')])
        printmsg('Level 1 job probably already estimated, moving on...', 'msgtitle', subnam{s}); continue;
    end
    printmsg(sprintf('Building Level 1 Job for %d Runs', length(rundir)),'msgtitle', subnam{s});

    % | Behavioral and Nuisance Regressor Files
    % | ========================================================================
    if ~isempty(nuisancepat)
        nuisance    = get_files([subdir{s} filesep runpat filesep nuisancepat], 'Nuisancer regressor')
    else
        nuisance = [];
    end
    behav       = get_files([subdir{s} filesep behavpat], 'Behavioral data');

    % | Get Images
    % | ========================================================================
    images          = cell(size(rundir));
    for r = 1:length(rundir)
        images{r} = files([rundir{r} filesep epipat]);
        if isempty(images{r})
            error('\nImage data not found! Failed search pattern:\n%s', [rundir{r} filesep epipat]);
        end
    end

    % | Run Loop
    % | ========================================================================
    for r = 1:length(rundir)

        % | Data for Current Run
        % | =====================================================================
        b = get_behavior(behav{r}, model, yesnokeys);
        b.blockwise(:,3) = b.blockwise(:,3) - adjons;

        % | Columns for b.blockwise
        % | =====================================================================
        % 1 - Block
        % 2 - Cond
        % 3 - Onset
        % 4 - Duration
        % 5 - Total_Errors
        % 6 - Foil_Errors

        % | Conditions
        % | =====================================================================
        for c = 1:length(b.condlabels)
            runs(r).conditions(c).name      = b.condlabels{c};
            runs(r).conditions(c).onsets    = b.blockwise(b.blockwise(:,2)==c, 3);
            runs(r).conditions(c).durations = b.blockwise(b.blockwise(:,2)==c, 4);
        end

        % | Floating Parametric Modulators
        % | =====================================================================
        if ~isempty(covidx)
            allpm           = b.blockwise(:,4:6);
            modelpm         = allpm(:,covidx);
            modelpmnames    = pmnames;
            novaridx = find(nanstd(modelpm)==0);
            if ~isempty(novaridx), modelpm(:,novaridx) = []; modelpmnames(novaridx) = []; end
            for p = 1:length(modelpmnames)
                runs(r).floatingpm(p).name = modelpmnames{p};
                runs(r).floatingpm(p).onsets = b.blockwise(:,3);
                runs(r).floatingpm(p).durations = b.blockwise(:,4);
                runs(r).floatingpm(p).values = modelpm(:,p);
            end
        end

    end


    if length(rundir)==1
        images = images{1};
        if iscell(nuisance), nuisance = nuisance{1}; end
    end


    % | General Information
    % | ========================================================================
    general_info.analysis           = analysisdir;
    general_info.is4D               = is4D;
    general_info.TR                 = TR;
    general_info.hpf                = HPF;
    general_info.autocorrelation    = armethod;
    general_info.nuisance_file      = nuisance;
    general_info.brainmask          = brainmask;
    general_info.maskthresh         = maskthresh;
    general_info.hrf_derivs         = [0 0];
    general_info.mt_res             = 16;
    general_info.mt_onset           = 8;

    % | Contrasts
    % | ========================================================================
    ncond   = length(b.condlabels);
    w1      = eye(ncond);
    if ncond==2
        w2 = [1 -1];
    else
        w2 = [.5 .5 -.5 -.5; .5 -5 .5 -.5; 1 0 -1 0; 0 1 0 -1; 1 -1 0 0; 0 0 1 -1; .5 -.5 -.5 .5];
    end
    weights = [w1; w2];
    ncon    = size(weights,1);
    for c = 1:ncon
        contrasts(c).type       = 'T';
        contrasts(c).weights    = weights(c,:);
        contrasts(c).name       = bspm_conweights2names(weights(c,:), b.condlabels);
    end
    if fcontrast
        contrasts(ncon+1).type      = 'F';
        contrasts(ncon+1).name      = 'Omnibus';
        contrasts(ncon+1).weights   = eye(ncond);
    end

    % | Make Job
    % | ========================================================================
    allinput{s} = bspm_level1(images, general_info, runs, contrasts);

    % | Cleanup Workspace
    % | ========================================================================
    clear general_info runs contrasts b modelpm modelpmnames

end

% | Run if desired
if run_it_now, bspm_runbatch(allinput); end

end
% =========================================================================
% * SUBFUNCTIONS
% =========================================================================
function b = get_behavior(in, opt, yesnokeys)
    % GET_BEHAVIOR
    %
    %   USAGE: b = get_behavior(in, opt)
    %
    %       in      behavioral data filename (.mat)
    %       opt     '2x2'  - full design
    %               '1x2'  - why vs. how
    %
    %       Columns for b.blockwise
    %          1 - Block
    %          2 - Cond
    %          3 - Onset
    %          4 - Duration
    %          5 - Total_Errors
    %          6 - Foil_Errors
    %
    % CREATED: Bob Spunt, Ph.D. (bobspunt@gmail.com) - 2014.02.24
    % =========================================================================
    if nargin < 1, error('USAGE: b = get_behavior(in, opt, yesnokeys)'); end
    if nargin < 2, opt = '2x2'; end
    if nargin < 3, yesnokeys = [1 2]; end
    if iscell(in), in = char(in); end

    % | read data
    % | ========================================================================
    d = load(in);
    % b.subjectID = d.subjectID;
    if ismember({'result'},fieldnames(d))
        data        = d.result.trialSeeker;
        blockwise   = d.result.blockSeeker;
    else
        data        = d.trialSeeker;
        blockwise   = d.blockSeeker;
    end

    % | blockwise accuracy and durations
    % | ========================================================================
    ntrials         = length(data(data(:,1)==1,1));
    data(data(:,8)==yesnokeys(1), 8) = 1;
    data(data(:,8)==yesnokeys(2), 8) = 2;
    data(:,10)      = data(:,4)~=data(:,8); % errors
    data(data(:,8)==0, 7:8) = NaN; % NR to NaN
    blockwise(:,3)  = data(data(:,2)==1, 6);
    blockwise(:,4)  = data(data(:,2)==ntrials, 9) - data(data(:,2)==1, 6);

    % | compute block-wise error counts
    % | ========================================================================
    for i = 1:size(blockwise, 1)
        blockwise(i,5) = sum(data(data(:,1)==i,10));  % all errors
        blockwise(i,6) = sum(data(data(:,1)==i & data(:,4)==2, 10)); % foil errors
    end

    % | re-code data
    % | ========================================================================
    con     = blockwise(:,2);
    condata = data(:,3);
    switch lower(opt)
        case {'2x2'}
            b.condlabels = {'Why_Face' 'Why_Hand' 'How_Face' 'How_Hand'};
        case {'1x2'}
            b.condlabels = {'Why' 'How'};
            blockwise(con==2, 2)                    = 1;
            data(condata==2, 3)                     = 1;
            blockwise(ismember(con,[3 4]), 2)       = 2;
            data(ismember(condata,[3 4]), 3)        = 2;
    end
    for i = 1:length(unique(data(:,3)))
        cdata = data(data(:,3)==i, [7 10]);
        b.accuracy(i)  = 100*(sum(cdata(:,2)==0)/size(cdata,1));
        b.rt(i)        = nanmean(cdata(:,1));
    end
    b.blockwise = blockwise;
    b.varlabels = {'Block' 'Cond' 'Onset' 'Duration' 'Total_Errors' 'Foil_Errors'};
end
function mfile_showhelp(varargin)
    % MFILE_SHOWHELP
    ST = dbstack('-completenames');
    if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
    eval(sprintf('help %s', ST(2).file));
end
function f = get_files(pat, label)
    f = files(pat);
    if isempty(f)
        error('\n\n%s file not found using search pattern:\n%s\n\n', label, pat)
    end
    fprintf('\n%s file:\n%s\n\n', label, char(f))
end
