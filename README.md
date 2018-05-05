# whyhowlocalizer-analysis
Software for analyzing single-subject fMRI data on the Why/How Localizer using SPM12

## Before you begin

- The main function to run on your end is `wrapper_level1_whyhow.m`. See below to learn about the various arguments you can pass to it to customize the analysis, and see the included `example_usage.m` file to run a model on real data (in the "example_data" folder).
- Make sure you have SPM12 and included subfolder `support` on your MATLAB search path before running
- If you plan to use Weighted Least Squares (WLS) estimation (I do this by default), you need to have a copy of the [Robust Weighted Least Squares toolbox](http://www.icn.ucl.ac.uk/motorcontrol/imaging/robustWLS.html) in the toolbox folder within your SPM12 directory. A copy of the toolbox is included in the file `rWLS_v4.0_SPM12.zip`.

## These varargins are used for finding the necessary files and directories

| NAME           | DEFINITION                                                                                 |
| :------------- | ------------------------------------------------------------------------------------------ |
| studydir       | full path to directory containing subject data folders                                     |
| behavpat       | search pattern for finding behavioral data within subject dirs, e.g., 'behav/whyhow*mat'   |
| epipat         | search pattern for finding functional data file(s) within each run, e.g., 'sw*nii'         |
| nuisancepat    | search pattern for finding nuisance regressor file within each run, e.g., 'rp*txt'         |
| runpat         | search pattern for finding run directories within subject dirs, e.g.,'raw/BOLD_WhyHow*'    |
| subpat         | search pattern for finding subject directories within studydir, e.g., 'Subject*'           |
| brainmask      | full path to brain mask to use (leave empty for none)                                      |

## These varargins are used to specify relevant details about the image and behavioral data being modeled

| NAME           | DEFINITION                                                                                 |
| :------------- | ------------------------------------------------------------------------------------------ |
| is4D           | flag for 4D image file (0=No, 1=Yes)                                                       |
| nskip          | number of initial TRs that have been removed (for adjusting stimulus onsets)               |
| TR             | acquisition repetition time (in seconds)                                                   |
| yesnokeys      | keys corresponding to yes/no responses (e.g., [1 2])                                       |

## These varargins are used to configure the model and estimation methods

| NAME           | DEFINITION                                                                                              |
| :------------- | ------------------------------------------------------------------------------------------------------- |
| basename       | base name for the analysis (e.g., 'WhyHow_SmoothedData')                                                |
| model          | string specifying model to use ('2x2' for full design, '1x2' to collapse faces/hands factor             |
| incl_err       | flag to include parametric covariate modeling blockwise variation in # of errors  (0=No, 1=Yes)         |
| incl_rt        | flag to include parametric covariate modeling blockwise variation in response time  (0=No, 1=Yes)       |
| armethod       | autocorrelation removal method (0=None, 1=AR(1), 2=Weighted Least Squares (WLS), 3=FAST)                |
| HPF            | high-pass filter cutoff to use (in seconds)                                                             |
| maskthresh     | implicit masking threshold (proportion of globals), default = 0.8                                       |
| fcontrast      | flag to compute omnibus F-contrast (useful for feature selection, e.g., in PPI analysis)  (0=No, 1=Yes) |
| run_it_now     | flag to run analysis now (0=No, 1=Yes)                                                                  |

