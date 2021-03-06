% Example Usage with a bunch of varargin

% add folders to MATLAB path
addpath(pwd); 
addpath('support');

% These varargins are used for finding the necessary files and directories
studydir     = fullfile(pwd, 'example_data'); 
subpat       = 'subject*';
behavpat     = 'behav/whyhow*mat';
runpat       = 'raw/BOLD_*';
epipat       = 'sw*nii';
nuisancepat  = 'rp*txt';

% These varargins are used to specify relevant details about the image and behavioral data being modeled
is4D        = 1;
nskip       = 4;

% These varargins are used to configure the model and estimation methods
run_it_now  = 1; 
incl_err    = 1;
incl_rt     = 1;
armethod    = 2;
HPF         = 100;

% Use wrapper to build the batch jobs
batchjob = wrapper_level1_whyhow( ...
        'run_it_now'  , run_it_now, ...
        'studydir'    , studydir, ...
        'subpat'      , subpat, ...
        'behavpat'    , behavpat, ...
        'runpat'      , runpat, ...
        'epipat'      , epipat, ...
        'nuisancepat' , nuisancepat, ...
        'is4D'        , is4D, ...
        'nskip'       , nskip, ...
        'incl_err'    , incl_err, ...
        'incl_rt'     , incl_rt, ...
        'armethod'    , armethod, ...
        'HPF'         , HPF);

