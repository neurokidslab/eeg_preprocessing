% robust_locdetrend() - Detrends the data in an EEGLAB EEG variable using
%                       local linear regression and a linear regression.
%                       Regression can be performed using ordinary least
%                       squares regression (which is fast) or a  method
%                       robust to outliers ("least angles regression"
%                       as implementd by the MATLAB function fit.m).
%                       The way local linear regression works is that a line is
%                       fit to all the data in a window, the window is then
%                       moved a step, and the process is repeated.  The
%                       lines from all windows are then averaged together
%                       via a weighted average to produce the estimated
%                       slow drift that is subtracted from the recording.
%
% Required Input:
%   EEG           - An EEGLAB EEG variable containing epoched data.
%
% Optional Inputs:
%   robust        - ['yes' or 'no'] If 'yes', least angles regression (a
%                   regression algorithm robust to outliers is used).
%                   Otherwise, ordinary least squares regression is used.
%                   {default: 'yes'}
%   wind_duration - The duration of the regression window in MILLISECONDS.
%                   For local linear regression, this should be less than
%                   the length of the epoch. The smaller this is, the
%                   higher the frequencies that will be affected by the
%                   filter. For strictly linear regression (i.e., NOT local
%                   linear regression), leave this empty or make it equal
%                   to the epoch length. {default: 0}
%   step_size     - The amount to increment the regression window at each
%                   step of the local linear regression process in
%                   MILLISECONDS.  The smaller this is, the higher the
%                   frequencies that will be affected by the filter.  As a
%                   rule of thumb, this should not exceed wind_duration/2.
%                   For strictly linear regression (i.e., NOT local linear
%                   regression), leave this empty or make it 0. {default:
%                   0)
%   make_plots    - ['yes' or 'no'] If 'yes', a plot will be made for each
%                   epoch and channel to illustrate the trend estimate
%                   and effect of detrending. The plot for an epoch will be
%                   held for two seconds and then it will be erased and
%                   replaced by the plot for the next epoch. {default:
%                   'no'}
%
% Outputs:
%   EEG           - Same as input EEG variable but with the following
%                   fields modified:
%                      EEG.data
%                      EEG.history
%                      EEG.saved
%                      EEG.icaact
%                      EEG.icspectra
%                      EEG.icfreqs
%
% Examples:
% -To perform robust linear detrending
% >>EEG=robust_locdetrend(EEG);
%
% -To perform ordinary least squares linear detrending
% >>EEG=robust_locdetrend(EEG,'no');
%
% -To perform robust local linear detrending with a one second window and
% a half second step between windows
% >>EEG=robust_locdetrend(EEG,'yes',1000,500);
%
% Notes:
%   When the length of the data is such that some data points at the end of
% the epoch don't fall in a regression window.  The last regression line is
% simply extended to the end of the data.
%
% Author:
% David Groppe
% Kutaslab, 3/2011

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 4/22/2011-'robust' and 'make_plots' options added
%

%%%%%%%%%%%%%%%% Future Work %%%%%%%%%%%%%%%%%
% Make work on continuous data
% Add verblevel?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ana Flo
% January 2018
% - It can be applied on non epoched data
% - If there is EEG.artifacts.BCT indicating where there are artifacts that
% data is ingnore

function EEG = robust_locdetrend(EEG,robust,wind_duration,step_size,make_plots)


%% Default Parameters
if nargin<2,
    robust=1; %1=use robust regression
elseif strcmpi(robust,'yes'),
    robust=1;
elseif strcmpi(robust,'no'),
    robust=0;
else
    error('Argument ''robust'' should be ''yes'' or ''no''.');
end
if (nargin<3)
    wind_duration=0;
end
if (nargin<4)
    step_size=0;
end
if (nargin<5)
    make_plots=0; %set to 1 to illustrate trend estimation for each epoch and channel of data
elseif strcmpi(make_plots,'yes'),
    make_plots=1;
elseif strcmpi(make_plots,'no'),
    make_plots=0;
else
    error('Argument ''make_plots'' should be ''yes'' or ''no''.');
end

%% Check Parameters
%conver wind_duration and step_size to seconds
wind_duration=wind_duration/1000;
step_size=step_size/1000;
if wind_duration<0,
    error('wind_duration needs to be greater than or equal to 0.');
end
if (step_size<0)
    error('step_size needs to be greater than or equal 0.');
elseif (step_size>wind_duration),
    error('step_size needs to be smaller than or equal to wind_duration.');
end


%% Update history of EEG variable
if ischar(EEG.history),
    temp=EEG.history;
    EEG.history=[];
    EEG.history{1}=temp;
end
EEG.history{length(EEG.history)+1}=['EEG=robust_locdetrend(EEG,' num2str(wind_duration) ...
    ',' num2str(step_size) ');'];

if ~wind_duration
    %doing strictly linear detrending (not local linear detrending)
    wind_duration=length(EEG.times)/EEG.srate;
    step_size=wind_duration;
end

% Are data continuous or epoched?
s=size(EEG.data);
if length(s)==3,
    continuous_data=0;
    nt=s(2); %number of time points per epoch
    Fs=EEG.srate;
    data_duration=nt./Fs;
    if (wind_duration>data_duration)
        error('wind_duration needs to be equal to or smaller than the duration of your epochs, %g.', ...
            data_duration);
    end
else
%     continuous_data=1;
    
    continuous_data=0;
    nt=s(2); %number of time points per epoch
    Fs=EEG.srate;
    data_duration=nt./Fs;
    if (wind_duration>data_duration)
        error('wind_duration needs to be equal to or smaller than the duration of your epochs, %g.', ...
            data_duration);
    end
    s(3) = 1;
end

if make_plots
    h1=figure;
end

%turn off robust regression warnings
warning('off','curvefit:fit:iterationLimitReached');
warning('off','curvefit:fit:nonDoubleYData');
if continuous_data,
    error('This function can''t handle continuous data yet.');
else
    n=wind_duration*Fs; % number of time points in window
    dn=step_size*Fs; % number of time points to increment window at each step
    xwt=((1:n)-n/2)/(n/2);
    wt=(1-abs(xwt).^3).^3; %weighting function for combining regression lines
    nwin=1+floor((nt-n)/dn); % number of time windows per epoch
    if ~nwin
        error('No time windows span these data.');
    end
    
    if robust,
        fo_ = fitoptions('method','LinearLeastSquares','Robust','LAR');
        ft_ = fittype('poly1');
    end
    
    for channel=1:s(1),
        fprintf('Robust detrending channel: %s\n',EEG.chanlocs(channel).labels);
        for epoch=1:s(3),
            y=EEG.data(channel,:,epoch)';
            
            % -------------------------------------------------------------
            % Put as zero periods wih artifacts
            if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BCT')
                bad = EEG.artifacts.BCT(channel,:,epoch)';
                if any(bad)
                    ini = find([bad(1); bad(2:end)-bad(1:end-1)]==1);
                    fin = find([bad(1:end-1)-bad(2:end); bad(end)]==1);
                    if ini(1)~=1 && fin(1)~=length(y)
                        y(ini(1):fin(1)) = (y(ini(1)-1)+y(fin(1)+1))/2;
                    elseif fin(1)~=length(y)
                        y(ini(1):fin(1)) = y(fin(1)+1);
                    else
                        y(ini(1):fin(1)) = 0;
                    end
                    if length(ini)>2
                        for iseg=2:length(ini)-1
                            y(ini(iseg):fin(iseg)) = (y(ini(iseg)-1)+y(fin(iseg)+1))/2;
                        end
                    end
                    if length(ini)>1
                        if fin(end)~=length(y)
                            y(ini(end):fin(end)) = (y(ini(end)-1)+y(fin(end)+1))/2;
                        elseif ini(end)~=1
                            y(ini(end):fin(end)) = y(ini(end)-1);
                        else
                            y(ini(end):fin(end)) = 0;
                        end
                    end
                end
%                 y(EEG.artifacts.BCT(channel,:,epoch)') = 0;
            end
            % -------------------------------------------------------------
            
            y_line=zeros(nt,1); % trend estimate
            norm=y_line; %sum of weights at each time point (we'll divide by this to average the linear fits in each window
            
            yfit=zeros(nwin,n); %linear fit for each moving window
            
            if make_plots
                figure(h1); clf;
                subplot(1,2,1);
                %plot raw data in Subplot 1
                plot([EEG.times(1) EEG.times(end)],[0 0],'k'); hold on;
                plot(EEG.times,y,'b');
            end
            
            for j=1:nwin
                use_time_pt_ids=dn*(j-1)+1:dn*(j-1)+n; %data within moving time window
                tseg=y(use_time_pt_ids);
                
                if robust,
                    cf_=fit([1:n]',tseg,ft_,fo_); %robust regression
                    a=cf_.p1; %slope
                    b=cf_.p2; %intercept
                else
                    %OLS regression parameters
                    y1=mean(tseg); %the mean of the segment of data in the window
                    y2=mean((1:n)'.*tseg)*2/(n+1);
                    a=(y2-y1)*6/(n-1); %a=slope,
                    b=y1-a*(n+1)/2;  %b=intercept
                end
                
                yfit(j,:)=(1:n)*a+b; %regression line
                y_line((j-1)*dn+(1:n))=y_line((j-1)*dn+(1:n))+(yfit(j,:).*wt)'; %weighted sum of regression lines
                norm((j-1)*dn+(1:n))=norm((j-1)*dn+(1:n))+wt'; %sum of the weights for combining regression lines
                
                if make_plots
                    figure(h1);
                    subplot(1,2,1);
                    hh=plot(EEG.times(use_time_pt_ids),yfit(j,:),'r');
                    set(hh,'linewidth',2);
                end
                
            end
            mask=find(norm>0);  %data points at the end might have a weight of 0, ignore them when normalizing
            y_line(mask)=y_line(mask)./norm(mask);
            indx=(nwin-1)*dn+n-1; %index of time points at the end that have a weight of 0
            npts=length(y)-indx+1;
            y_line(indx:end)=(n+1:n+npts)'*a+b; %flesh out final points with last estimated slope & intercept
            
            EEG.data(channel,:,epoch)=[EEG.data(channel,:,epoch)'-y_line]'; %detrended data
%             EEG.data(channel,:,epoch)=[y-y_line]'; %detrended data
            
            if make_plots
                figure(h1);
                subplot(1,2,1);
                hh=plot(EEG.times,y_line,'g--');
                set(hh,'linewidth',2);
                v=axis();
                axis([EEG.times(1) EEG.times(end) v(3:4)]);
                [LEGH, OBJH]=legend('EEG','Local Linear Fits','Weighted Mean of Fits','location','best');
                set(OBJH(4),'color','b','linewidth',2);
                set(OBJH(6),'color','r')
                set(OBJH(8),'color','g','linestyle','--');
                hh=xlabel('Time (ms)');
                hh=ylabel('\muV');
                title('Raw Data');
                
                %plot detrended data in Subplot 2
                subplot(1,2,2);
                plot([EEG.times(1) EEG.times(end)],[0 0],'k'); hold on;
                hh=plot(EEG.times,EEG.data(channel,:,epoch),'b'); %detrended data
                %v=axis();
                axis([EEG.times(1) EEG.times(end) v(3:4)-mean(v(3:4))]);
                hh=xlabel('Time (ms)');
                hh=ylabel('\muV');
                title('Local Linear Detrended Data');
                
                pause(2);
            end
        end
    end
end
%turn robust regression warnings back on
warning('on','curvefit:fit:iterationLimitReached');
warning('on','curvefit:fit:nonDoubleYData');

EEG.saved='no';

% Recompute ICA activations and remove IC spectral info if they exist
if ~isempty(EEG.icaact),
    EEG=recomp_act(EEG);
    
    if isfield(EEG,'icspectra'),
        EEG.icspectra=[];
    end
    
    if isfield(EEG,'icfreqs'),
        EEG.icfreqs=[];
    end
end

function EEG=recomp_act(EEG)
%function EEG=recomp_act(EEG)
%
% Recompute ICA activations

s=size(EEG.data);
if length(s)==3,
    %epoched data
    EEG.icaact=reshape(EEG.data,s(1),s(2)*s(3));
    EEG.icaact=EEG.icaweights*EEG.icasphere*EEG.icaact;
    EEG.icaact=reshape(EEG.icaact,s(1),s(2),s(3));
else
    %continuous data
    EEG.icaact=EEG.icaweights*EEG.icasphere*EEG.data;
end


