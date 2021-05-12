% Function that identifies motion artifacts. by appling different
% algorithms.
% The steps to be applied for the identification are indicated in the 
% structure Art
%
% -------------------------------------------------------------------------

function EEG = eega_tArtifacts( EEG, Art, varargin )

% possible steps to detect artifacts
StepsDetect = { 'eega_tRejPwr',...
                'eega_tRejCorrCh',...
                'eega_tRejTimeVar',...
                'eega_tRejAmp',...
                'eega_tRejRunningAvg',...
                'eega_tRejFastChange',...
                'eega_tRejDerivate',...
                'eega_tRejAmpElecVar',...
                };
% possible steps that are appplied always at the end of each loop
StepsPostDetect = { 'eega_tRejChPercSmpl';...
                    'eega_tRejSmplPercCh';...
                    'eega_tIncShortBad';...
                    'eega_tRejShortGood';...
                    'eega_tMask'};

%% ========================================================================
%% Parameters

% default parameters
cnf.FilterLp        = [];
cnf.FilterHp        = [];
cnf.KeepRejPre      = 1;
cnf.KeepRejCause    = 0;
cnf.MaxLoop         = [];
cnf.RejTolerance    = 0;
cnf.Lim2RejSbj      = 1;

% get extra inputs
[cnf, OK, extrainput] = eega_getoptions(cnf, varargin);
if ~OK
    error('eega_tArtifacts: Non recognized inputs')
end

% check the Art input structure
Steps = cell(length(Art),1);
ArtLoops = cell(length(Art),1);
for i=1:length(Art)
    if isfield(Art(i),'algorithm') && ~isempty(Art(i).algorithm)
        Steps{i} = Art(i).algorithm;
    else
        error('eega_tArtifacts: Field ".algorithm" is missing for step %d.',i)
    end
    if isfield(Art,'loops')
        ArtLoops{i} = Art(i).loops;
    else
        error('eega_tArtifacts: Field ".loop" is missing for step %d.',i)
    end
    if ~isfield(Art(i),'P')
        error('eega_tArtifacts: Field ".P" is missing for step %d.',i)
    end
end

% determine the maximun number of loops
if isempty(cnf.MaxLoop)
    theloops = [];
    for i=1:length(ArtLoops)
        theloops = [theloops(:); ArtLoops{i}(:)];
    end
    cnf.MaxLoop = unique(max(theloops));
end

cnf.FilterDo = ~isempty(cnf.FilterLp) | ~isempty(cnf.FilterHp);

%% ========================================================================
%% Show the info reharding the parameters
nEl = size(EEG.data,1);
nS = size(EEG.data,2);
nEp = size(EEG.data,3);
nsmpl = nEl*nS*nEp;
nSteps = length(Art);

fprintf('\n')
fprintf('### MOTION ARTIFACT REJECTION ALGORITHMS ###\n');
fprintf('\n')
fprintf('%d electrodes, %d samples, %d epochs\n', nEl, nS, nEp);
fprintf('Initial number of samples: %d\n', nEl*nS*nEp);
fprintf('\n')
if cnf.KeepRejCause==0
    fprintf('- The rejection cause will not be saved \n')
else
    fprintf('- The rejection cause will be saved \n')
end
if cnf.KeepRejPre==0
    fprintf('- BCT will be reset \n')
else
    fprintf('- The previous BCT will be kept \n')
end
if ~isempty(cnf.FilterLp)
    fprintf('- Data will be low-pass filtered before detection: %5.2f Hz\n',cnf.FilterLp)
end
if ~isempty(cnf.FilterHp)
    fprintf('- Data will be low-pass filtered before detection: %5.2f Hz\n',cnf.FilterHp)
end
fprintf('- The rejection algorithms will be applied: \n')
fprintf('   --> a maximun of %d times \n', cnf.MaxLoop)
fprintf('   --> or till the new rejected data is less than %6.2f %% \n', cnf.RejTolerance)
fprintf('\n')

%% ========================================================================
%% Get already rejected data if the rejection structure already exist
if ~isfield(EEG, 'artifacts')
    EEG.artifacts.algorithm.parameters = [];
    EEG.artifacts.algorithm.stepname = [];
    EEG.artifacts.algorithm.rejxstep = [];
    EEG.artifacts.BCT = false(nEl,nS,nEp);
    EEG.artifacts.BC = false(nEl,1,nEp);
    EEG.artifacts.BCmanual = [];
    EEG.artifacts.BT = false(1,nS,nEp);
    EEG.artifacts.BE = false(1,1,nEp);
    EEG.artifacts.BS = false(1);
end
if cnf.KeepRejPre==0
    EEG.artifacts.algorithm.parameters = [];
    EEG.artifacts.algorithm.stepname = [];
    EEG.artifacts.algorithm.rejxstep = [];
    EEG.artifacts.BCT = false(nEl,nS,nEp);
    EEG.artifacts.BCT = false(nEl,nS,nEp);
    EEG.artifacts.BC = false(nEl,1,nEp);
    EEG.artifacts.BT = false(1,nS,nEp);
    EEG.artifacts.BE = false(1,1,nEp);
    EEG.artifacts.BS = false(1);
end

BCT = EEG.artifacts.BCT;
BCTr = cell(nSteps,1);
BCTr{1} = BCT;

idxpost = false(1,length(Art));
for i=1:length(Art)
    iii = strcmp(Art(i).algorithm,StepsPostDetect(:));
    if any(iii)
        idxpost(i)=1;
    end
end
if any(idxpost)
    StepsPost = Steps(idxpost);
else
    StepsPost = {};
end
if any(~idxpost)
    StepsArt = Steps(~idxpost);
else
    StepsArt = {};
end

%% ========================================================================
%% Filter data before rejection if requeired
if cnf.FilterDo    
    dat_orig = EEG.data;
    if ~isempty(cnf.FilterLp)
        EEG = pop_eegfiltnew(EEG, [], cnf.FilterLp, [], 0, [], [], 0);
    end
    if ~isempty(cnf.FilterHp)
        EEG = pop_eegfiltnew(EEG, cnf.FilterHp, [], [], 0, [], [], 0);
    end
end

%% ========================================================================
%% Reject based on the different algorithms
stepname = cell(nSteps,1);
stepdone = false(nSteps,1);
parameters = cell(nSteps,1);
RejxStep = zeros(nSteps,1);
RejxStepNew = zeros(nSteps,1);
ok = 0;
loop = 0;

while ~ok
    
    loop = loop+1;    
    fprintf('\n******** LOOP %d ********\n\n',loop)
    BCT_pre = BCT;
        
    %% --------------------------------------------------------------------
    %% Algorithms to detect artifacts
    for i=1:length(Art) 
        if any(strcmp(StepsArt,Art(i).algorithm)) && any(Art(i).loops==loop)
            step = sum(stepdone)+1;
            thealgh = Art(i).algorithm;
            flds = fieldnames(Art(i).P);
            vals = struct2cell(Art(i).P);
            inputs = cat(1,flds(:)',vals(:)');
            fhandle = str2func(thealgh);
            
            fprintf('%s\n',repmat('-',[1 30]))
            fprintf('Rejection step % 2d: %s\n\n',step,thealgh)
            
            [ ~, bct ] = fhandle( EEG, inputs{:});
            
            
            nsmplrej = sum(bct(:));
            nnewsmplrej = bct & ~BCT;
            nnewsmplrej = sum(nnewsmplrej(:));
            RejxStep(step) = nsmplrej;
            RejxStepNew(step) = nnewsmplrej;
            BCT = bct | BCT;
            BCTr{step} = bct;
            
            fprintf('New data rejected % 4.2f %%\n\n', nnewsmplrej/nsmpl*100)
            
            stepdone(step) = 1;
            stepname{step} = thealgh;
            parameters{step} = Art(i).P;
            
        end
    end
    
    % Update rejection
    EEG.artifacts.BCT = BCT;
    
    
    %% --------------------------------------------------------------------
    %% Algorithm to reject more data based on the rejected data
    for i=1:length(Art) 
        if any(strcmp(StepsPost,Art(i).algorithm)) && any(Art(i).loops==loop)
            step = sum(stepdone)+1;
            thealgh = Art(i).algorithm;
            flds = fieldnames(Art(i).P);
            vals = struct2cell(Art(i).P);
            inputs = cat(1,flds(:)',vals(:)');
            fhandle = str2func(thealgh);
            
            fprintf('%s\n',repmat('-',[1 30]))
            fprintf('Rejection step % 2d: %s\n\n',step,thealgh)
            
            [ EEG, change ] = fhandle( EEG, inputs{:},'updatesummary',0,'updatealgorithm',0);
            
            nsmplchange = sum(change(:));
            RejxStep(step) = nsmplchange;
            RejxStepNew(step) = nsmplchange;
            BCTr{step} = change;
            if strcmp(thealgh,'eega_tIncShortBad')
                fprintf('New data re-included % 4.2f %%\n\n', nsmplchange/nsmpl*100)
            else
                fprintf('New data rejected % 4.2f %%\n\n', nsmplchange/nsmpl*100)          
            end
            
            stepdone(step) = 1;
            stepname{step} = thealgh;
            parameters{step} = Art(i).P;
        end
    end
    
    % Update rejection
    BCT = EEG.artifacts.BCT;
    
  
    
    %% --------------------------------------------------------------------
    %% See if new data was rejected in this loop
    newrej = BCT & ~BCT_pre;
    newrej = sum(newrej(:))/(nEl*nS*nEp)*100;
    fprintf('%s\n',repmat('-',[1 30]))
    fprintf('NEW DATA REJECTED LOOP %d: %8.4f %%\n',loop, newrej)
    fprintf('%s\n',repmat('-',[1 30]))
    if (cnf.RejTolerance ~= 0) && (newrej <= cnf.RejTolerance)
        ok=1;        
    end
    if loop==cnf.MaxLoop
        ok=1;
    end
      
    
end

fprintf('\nEnd of loops \n')

%% ========================================================================
%% Original data back
if cnf.FilterDo
    EEG.data = dat_orig;
end

%% ========================================================================
%% Summary 
step_tot = sum(stepdone);
TotSmpl = nEl*nS*nEp;
TotSmplRej = sum(BCT(:));
TotSmplRem = TotSmpl-TotSmplRej;
ProRejxStep = RejxStep / TotSmpl*100;
ProRejxStepNew = RejxStepNew / TotSmpl*100;
fprintf('%s\n',repmat('\',[1 80]))
for i=1:step_tot
    stp = stepname{i};
    fprintf('  - Step %02d: %s %s %12d (% 5.1f %%) / New %12d (% 5.1f %%) \n',...
        i, stp, repmat(' ',[20-numel(stp) 1]), RejxStep(i), ProRejxStep(i), RejxStepNew(i), ProRejxStepNew(i));
end
fprintf('%s\n',repmat('\',[1 80]))
fprintf('Rejected samples:  %010d (%2.1f %%)\n', TotSmplRej, TotSmplRej/TotSmpl*100);
fprintf('Remaining samples: %010d (%2.1f %%)\n', TotSmplRem, TotSmplRem/TotSmpl*100);
fprintf('%s\n',repmat('\',[1 80]))
fprintf('\n')

%% ========================================================================
%% Update EEG
EEG.artifacts.algorithm.parameters = cat(1,EEG.artifacts.algorithm.parameters(:),parameters(:));
stepname = stepname(logical(stepdone));
EEG.artifacts.algorithm.stepname = cat(1,EEG.artifacts.algorithm.stepname(:),stepname(:));
EEG.artifacts.algorithm.rejxstep = cat(1,EEG.artifacts.algorithm.rejxstep(:),RejxStep(:));
if exist('eega_summarypp','file')==2
    EEG = eega_summarypp(EEG);
end
if cnf.KeepRejCause
    EEG.artifacts.BCTSr = BCTr(logical(stepdone));
end

end
