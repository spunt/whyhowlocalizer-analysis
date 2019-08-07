function b = score_whyhow_socns(in, varargin)
% SCORE_WHYHOW_SOCNS  Score Yes/No Why/How Task Behavior
%
%  USAGE: b = score_whyhow_socns(in, varargin)
% __________________________________________________________________________
%  INPUT
%	in:             behavioral data MAT-file
%
% ________________________________________________________________________________________
%  VARARGIN (partial matches OK; run without arguments to see default values)
% | NAME            | DESCRIPTION                                                        |
% |-----------------|--------------------------------------------------------------------|
% | taskModel       | '2x2' for why/how:soc/ns, '2x3' for full design                   |
% | yesNoKeys       | keys subject pressed for yes/no                                    |
% | omitNoResponses | flag to omit no response trials                                    |
% | omitBadBlocks   | flag to omit blocks with errors on >= badBlockThresh trials        |
% | badBlockThresh  | threshold for omitBadBlocks                                        |
% | omitOutlierRTs  | flag to omit responses with outlier RTs defined by outlierSDThresh |
% | outlierSDThresh | SD threshold for                                                   |
% | omitTooFast     | flag to omit response with implausibly fast RTs                    |
% | tooFastThresh   | threshold for omitTooFast (seconds)                                |
% ________________________________________________________________________________________
% b.blockwise - column key
%     1 - Block #
%     2 - Cond: 1=Why_NS, 2=Why_Face, 3=Why_Hand, 4=How_NS, 5=How_Face, 6=How_Hand
%     3 - Block Onset (secs): corresponds to onset of first trial
%     4 - Block Duration (secs): corresponds to offset of last trial
%     5 - # Errors
%     6 - # Errors (Foils)
%     7 - # No Response
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-03-06
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
def  = { ...
  'taskModel'           ,       '2x3'                       ,...
  'yesNoKeys'           ,       [1 2]                       ,...
  'omitNoResponses'     ,       1                           ,...
  'omitBadBlocks'       ,       0                           ,...
  'badBlockThresh'      ,       4                           ,...
  'omitOutlierRTs'      ,       0                           ,...
  'outlierSDThresh'     ,       4                           ,...
  'omitTooFast'         ,       1                           ,...
  'tooFastThresh'       ,       .25                          ...
};
vals = setargs(def, varargin);
if nargin<1, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if iscell(in), in = char(in); end

% | CONFIGURATION STRING
b.defstr = sprintf('%s_Outlier%d_NR%d_BadBlocks%d_TooFast%d', taskModel, omitOutlierRTs, omitNoResponses, omitBadBlocks, omitTooFast);

% | read data
% | ========================================================================
d = load(in);
b.subjectID = d.subjectID;
if ismember({'result'},fieldnames(d))
    data        = d.result.trialSeeker;
    blockwise   = d.result.blockSeeker;
    qidx        = blockwise(:, end);
    questions   = regexprep(d.result.preblockcues(qidx), 'Is the person ', '');
else
    data        = d.trialSeeker;
    blockwise   = d.blockSeeker;
    questions   = d.ordered_questions;
end
strrm = {'Is the photo ' 'Is it a result of a ' 'Is it a result of ' 'Is it going to result in a ' 'Is the person '};
for i = 1:length(strrm)
    questions = regexprep(questions, strrm{i}, '');
end
questions = regexprep(questions, ' ', '_');
questions = regexprep(questions, '?', '');
questions = regexprep(questions, '-', '_');

% | blockwise accuracy and durations
% | ========================================================================
ntrials                         = length(data(data(:,1)==1,1));
actresp                         = data(:,8);
data(actresp==yesNoKeys(1), 8)  = 1;
data(actresp==yesNoKeys(2), 8)  = 2;
data(:,10)          = data(:,4)~=data(:,8); % error trials (includes no response here)
SKIPIDX             = data(:,8)==0;
data(SKIPIDX, 7:8)  = NaN; % converts No Response to NaN
if omitNoResponses, data(SKIPIDX,3) = NaN; end
if omitTooFast
    FASTIDX           = data(:,7) < tooFastThresh;
    data(FASTIDX,[3 7 8]) = NaN;
end
blockwise(:,3)      = data(data(:,2)==1, 6);
blockwise(:,4)      = data(data(:,2)==ntrials, 9) - data(data(:,2)==1, 6);
nskip               = sum(SKIPIDX);
nerror              = sum(data(:,10));
nfoilerror          = sum(data(data(:,4)==2, 10));
b.pctskip           = 100*(nskip/size(data,1));
b.pcterror          = 100*(nerror/size(data,1));
b.pctfoilerror      = 100*(nfoilerror/sum(data(:,4)==2));

% | compute block-wise error counts
% | ========================================================================
for i = 1:size(blockwise, 1)
    blockwise(i,5) = sum(data(data(:,1)==i,10));  % all errors
    blockwise(i,6) = sum(data(data(:,1)==i & data(:,4)==2, 10)); % foil errors
    blockwise(i,7) = sum(isnan(data(data(:,1)==i,7))); % no response
end
if omitBadBlocks
    BADBLOCK = blockwise(blockwise(:,end)>=badBlockThresh, 1);
    data(ismember(data(:,1), BADBLOCK), [3 7 8]) = NaN;
end

% | Trialwise Data
% | ========================================================================
desmatfile = fullfile(fileparts(mfilename('fullpath')), 'design_socns.mat');
if exist(desmatfile)
    des = load(fullfile(fileparts(mfilename('fullpath')), 'design_socns.mat'));
    qim = des.alldesign{1}.qim(data(:,5),:);
    b.trialwise = [qim num2cell(data(:,[3 4 7 10]))];
    b.trialwise = sortrows(b.trialwise, [2 1]);
end

% | re-code data
% | ========================================================================
con     = blockwise(:,2);
condata = data(:,3);
switch taskModel
    case {'2x3'}
        b.condlabels    = {'Why-NS' 'Why-Face' 'Why-Hand' 'How-NS' 'How-Face' 'How-Hand'};
        qcond = b.condlabels(blockwise(:,2))';
        rt = round(blockwise(:,4)*1000);
        err = blockwise(:,5);
        b.blocklabels = strcat(upper(qcond), '-', upper(questions), '-', num2str(rt), 'ms', '-', num2str(err), 'error');
        b.blocklabels = regexprep(b.blocklabels, ' ', '');
        wbnames         = {'WhyRTcost_NS' 'WhyRTcost_FACE' 'WhyRTcost_HAND'};
    case {'2x2'}
        b.condlabels    = {'Why_NS' 'Why_Soc' 'How_NS' 'How_Soc'};
        wbnames         = {'WhyRTcost_NS' 'WhyRTcost_SOC'};
        blockwise(ismember(con,[2 3]), 2)           = 2;
        blockwise(con==4, 2)                        = 3;
        blockwise(ismember(con,[5 6]), 2)           = 4;
        condata(ismember(condata,[2 3]))            = 2;
        condata(condata==4)                         = 3;
        condata(ismember(condata,[5 6]))            = 4;
end
data(:,3)   = condata;
data(isnan(data(:,3)), :) = [];
ncond       = length(unique(data(:,3)));
b.rtcell    = cell(ncond, 1);
for i = 1:ncond
    cdata           = data(data(:,3)==i, [7 10]);
    tmpacc          = data(data(:,3)==i, [4 8]);
    tmpacc(tmpacc==2) = 0;
    [b.dprime(i), b.responsebias(i)] = calc_signaldetectiontheory(tmpacc);
    b.acccell{i}    = tmpacc;
    b.accuracy(i)   = 100*(sum(cdata(:,2)==0)/size(cdata,1));
    tmprt           = cdata(cdata(:,2)==0, 1);
    if all([omitOutlierRTs ~isempty(tmprt)])
        tmprt = outlier2nan(tmprt, outlierSDThresh);
    end
    b.rtcell{i}     = tmprt;
    b.rt(i)         = nanmean(tmprt);
    if b.accuracy(i)==0, b.accuracy(i) = NaN; end
end
b.blockwise = blockwise;
b.varlabels = {'Block' 'Cond' 'Onset' 'Duration' 'Total_Errors' 'Foil_Errors'};
for i = 1:(ncond/2)
    rtpooledsd  = nanstd(vertcat(b.rtcell{[i i+(ncond/2)]}));
    whym = nanmean(b.rtcell{i});
    howm = nanmean(b.rtcell{i+(ncond/2)});
    b.whyrtcost(i) = (howm - whym)/rtpooledsd;
end

varnames    = {'ACC' 'RT' 'DPRIME' 'BIAS'};
% respcat     = {'All' 'Yes' 'No'};
condnames   = b.condlabels;
allnames    = [];

for i = 1:length(varnames)
    allnames = [allnames strcat(varnames{i}, '_', condnames)];
end

accnames    = [];
rtnames     = [];
for c = 1:length(condnames)
    tmp = {['ACC_' condnames{c}]};
    accnames = [accnames tmp];
    tmp = {['RT_' condnames{c}]};
    rtnames = [rtnames tmp];
end
b.allname = [allnames wbnames];
b.alldata = horzcat(b.accuracy, b.rt, b.dprime, b.responsebias, b.whyrtcost);
b = orderfields(b);
end
function [dprime, responsebias] = calc_signaldetectiontheory(tmpacc)

true_rescale        = 100/sum(tmpacc(:,1)==1);
false_rescale       = 100/sum(tmpacc(:,1)==0);
count.hits          = true_rescale*(sum(tmpacc(:,1)==1 & tmpacc(:,2)==1)) + 1;
count.misses        = true_rescale*sum(tmpacc(:,1)==1 & tmpacc(:,2)==0) + 1;
count.falsealarms   = false_rescale*sum(tmpacc(:,1)==0 & tmpacc(:,2)==1) + 1;
count.correctreject = false_rescale*sum(tmpacc(:,1)==0 & tmpacc(:,2)==0) + 1;
hit_rate            = count.hits / (count.hits + count.misses);
zhr                 = norminv(hit_rate);
false_alarm_rate    = count.falsealarms / (count.falsealarms + count.correctreject);
zfar                = norminv(false_alarm_rate);
% | D-Prime
dprime              = zhr - zfar;
% | Response Bias (Criterion, c)
responsebias           = -(zhr + zfar)/2;

end
function argstruct = setargs(defaultargs, varargs)
    % SETARGS Name/value parsing and assignment of varargin with default values
    %
    % This is a utility for setting the value of optional arguments to a
    % function. The first argument is required and should be a cell array of
    % "name, default value" pairs for all optional arguments. The second
    % argument is optional and should be a cell array of "name, custom value"
    % pairs for at least one of the optional arguments.
    %
    %   USAGE: argstruct = setargs(defaultargs, varargs)
    % __________________________________________________________________________
    % OUTPUT
    %
    % 	ARGSTRUCT
    %    structure containing the final argument values
    % __________________________________________________________________________
    % INPUTS
    %
    % 	DEFAULTARGS
    %     cell array of "'Name', value" pairs for all variables with default
    %     values
    %
    % 	VARARGS [optional]
    %     cell array of user-specified "'Name', value" pairs for one or more of
    %     the variables with default values. this will typically be the
    %     "varargin" cell array. for each pair, SETARGS determines if the
    %     specified variable name can be uniquely matched to one of the default
    %     variable names specified in DEFAULTARGS. matching uses STRNCMPI and
    %     thus is case-insensitive and open to partial name matches (e.g.,
    %     default variable name 'FontWeight' would be matched by 'fontweight',
    %     'Fontw', etc.). if a match is found, the user-specified value is then
    %     used in place of the default value. if no match is found or if
    %     multiple matches are found, SETARGS returns an error and displays in
    %     the command window information about the argument that caused the
    %     problem.
    % __________________________________________________________________________
    % USAGE EXAMPLE (TO BE USED AT TOP OF FUNCTION WITH VARARGIN)
    %
    %     defaultargs = {'arg1', 0, 'arg2', 'words', 'arg3', rand};
    %     argstruct   = setargs(defaultargs, varargin)
    %


    % ---------------------- Copyright (C) 2015 Bob Spunt -----------------------
    %	Created:  2015-03-11
    %	Email:    spunt@caltech.edu
    %
    %   This program is free software: you can redistribute it and/or modify
    %   it under the terms of the GNU General Public License as published by
    %   the Free Software Foundation, either version 3 of the License, or (at
    %   your option) any later version.
    %       This program is distributed in the hope that it will be useful, but
    %   WITHOUT ANY WARRANTY; without even the implied warranty of
    %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    %   General Public License for more details.
    %       You should have received a copy of the GNU General Public License
    %   along with this program.  If not, see: http://www.gnu.org/licenses/.
    %
    if nargin < 1, mfile_showhelp; return; end
    if nargin < 2, varargs = []; end
    defaultargs = reshape(defaultargs, 2, length(defaultargs)/2)';
    if ~isempty(varargs)
        if mod(length(varargs), 2)
            error('Optional inputs must be entered as "''Name'', Value" pairs, e.g., myfunction(''arg1'', val1, ''arg2'', val2)');
        end
        arg = reshape(varargs, 2, length(varargs)/2)';
        for i = 1:size(arg,1)
        idx = strncmpi(defaultargs(:,1), arg{i,1}, length(arg{i,1}));
        if sum(idx) > 1
            error(['Input "%s" matches multiple valid inputs:' repmat('  %s', 1, sum(idx))], arg{i,1}, defaultargs{idx, 1});
        elseif ~any(idx)
            error('Input "%s" does not match a valid input.', arg{i,1});
        else
            defaultargs{idx,2} = arg{i,2};
        end
        end
    end
    for i = 1:size(defaultargs,1), assignin('caller', defaultargs{i,1}, defaultargs{i,2}); end
    if nargout>0, argstruct = cell2struct(defaultargs(:,2), defaultargs(:,1)); end
    end
    function mfile_showhelp(varargin)
    % MFILE_SHOWHELP
    ST = dbstack('-completenames');
    if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
    eval(sprintf('help %s', ST(2).file));
    end
function y = outlier2nan(x, sd, rmcase)
    % OUTLIER2NAN
    %
    % Find outliers recusviely using leave-one-out method and convert to NaN.
    % If input X is a matrix, outliers computed columnwise.
    %
    %   USAGE: y = outlier2nan(x, sd, rmcase)
    % _________________________________________________________________________
    % INPUTS
    %   x       = input matrix
    %   sd      = number of standard deviations to truncate
    %   rmcase  = option to remove case/row
    %
        if nargin < 1, mfile_showhelp; return; end
        if nargin < 2, sd = 3; end
        if nargin < 3, rmcase = 0; end

        y       = zeros(size(x));
        ncol    = size(x,2);
        idx2y   = 1:ncol;

        while ~isempty(x)

            nnan1                      = sum(isnan(x));
            x(oneoutzscore(x, 1) > sd) = NaN;
            nancol                     = (sum(isnan(x)) - nnan1)==0;
            y(:, idx2y(nancol))        = x(:, nancol);
            x(:, nancol)               = [];
            idx2y(nancol)              = [];

        end

        if rmcase, y(isnan(mean(y,2)),:) = []; end

    end
function zx = oneoutzscore(x, returnas)
        % ONEOUTZSCORE Perform columnwise leave-one-out zscoring
        %
        % USAGE: zx = oneoutzscore(x, returnas)
        %
        %   returnas: 0, signed values (default); 1, absolute values
        %
        if nargin<1, disp('USAGE: zx = oneoutzscore(x, returnas)'); return; end
        if nargin<2, returnas = 0; end
        if size(x,1)==1, x=x'; end
        zx              = x;
        [nrow, ncol]    = size(x);
        for c = 1:ncol
            cin         = repmat(x(:,c), 1, nrow);
            theoneout   = cin(logical(eye(nrow)))';
            theleftin   = reshape(cin(logical(~eye(nrow))),nrow-1,nrow);
            cz          = (theoneout-nanmean(theleftin))./nanstd(theleftin);
            zx(:,c)     = cz';
        end
        if returnas, zx = abs(zx); end
    end
