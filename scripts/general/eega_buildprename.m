function [Fnct, fnctname] = eega_buildprename(Fnct, do)

%% remove functions

% all the functions
thef = Fnct;
thef(2:2:end) = [];
thef = thef(:);

% functions options
if ~isempty(do)
    fopt = fieldnames(do);
    s_rmv = zeros(length(Fnct),1);
    for i=1:length(fopt)
        fopti = fopt{i};
        if any(strcmp(thef,fopti)) && ~do.(fopti)
            % find the postition of each function
            idx = find(strcmp(thef,fopti));
            idx = [idx*2-1; idx*2];
            s_rmv(idx) = 1;
        end
    end
    Fnct(logical(s_rmv)) = [];
end

%% Determine the pre-name for the file
fnctname = [];
for i=1:2:length(Fnct)
    new='';
    switch Fnct{i}
        case 'eega_tTargetPCAxElEEG'
            new = 'pca_';
        case 'eega_tInterpSpatialEEG'
            new = 'i_';
        case 'eega_tInterpSpatialSegmentEEG'
            new = 'is_';
%             if any(strcmp(Fnct{i+1},'splicemethod')) && Fnct{i+1}{find(strcmp(Fnct{i+1},'splicemethod'))+1}
%                 new = sprintf('i%d_',Fnct{i+1}{find(strcmp(Fnct{i+1},'splicemethod'))+1});
%             else
%                 new = 'i0_';
%             end
        case 'eega_tArtifacts'
            new = 'a_';
        case 'eega_detrend'
            new = 'd_';
        case 'eega_ntdetrend'
            new = 'rd_';
        case 'pop_resample'
            new = sprintf('rs%d_',Fnct{i+1}{1});
        case 'eega_resample'
            new = sprintf('rs%d_',Fnct{i+1}{1});
        case 'eega_filter'
            if any(strcmp(Fnct{i+1},'minphase')) && Fnct{i+1}{find(strcmp(Fnct{i+1},'minphase'))+1}
                causal = 'c';
            else
                causal = '';
            end
            new = [];
            if ~isempty(Fnct{i+1}{1}) && isempty(Fnct{i+1}{2})
                new = [new sprintf('fh%s%1.0d_', causal, Fnct{i+1}{1})];
            end
            if ~isempty(Fnct{i+1}{2}) && isempty(Fnct{i+1}{1})
                new = [new sprintf('fl%s%d_', causal, Fnct{i+1}{2})];
            end
            if ~isempty(Fnct{i+1}{1}) && ~isempty(Fnct{i+1}{2})
                if length(Fnct{i+1})>2
                    idx = find(strcmp(Fnct{i+1}(3:end), 'revfilt'));
                    if ~isempty(idx) && Fnct{i+1}{2+idx+1}==1
                        notch = 1;
                    else
                        notch = 0;
                    end
                else
                    notch = 0;
                end
                if notch
                    new = [new sprintf('fn%s%d-%d_', causal, Fnct{i+1}{1}, Fnct{i+1}{2})];
                else
                    new = [new sprintf('fb%s%1.0d-%d_', causal, Fnct{i+1}{1}, Fnct{i+1}{2})];
                end
            end
        case 'pop_eegfiltnew'
            if ~isempty(Fnct{i+1}{7}) && Fnct{i+1}{7}==1  % minphase==1
                causal = 'c';
            else
                causal = '';
            end
            new = [];
            if ~isempty(Fnct{i+1}{1}) && isempty(Fnct{i+1}{2})
                new = [new sprintf('fh%s%1.0d_', causal, Fnct{i+1}{1})];
            end
            if ~isempty(Fnct{i+1}{2}) && isempty(Fnct{i+1}{1})
                new = [new sprintf('fl%s%d_', causal, Fnct{i+1}{2})];
            end
            if ~isempty(Fnct{i+1}{1}) && ~isempty(Fnct{i+1}{2})
                if ~isempty(Fnct{i+1}{4}) && Fnct{i+1}{4}  % revfilt==1
                    notch = 1;
                else
                    notch = 0;
                end
                if notch
                    new = [new sprintf('fn%s%d-%d_', causal, Fnct{i+1}{1}, Fnct{i+1}{2})];
                else
                    new = [new sprintf('fb%s%1.0d-%d_', causal, Fnct{i+1}{1}, Fnct{i+1}{2})];
                end
            end
        case 'eega_epoch'
            if length(Fnct{i+1}{1})==1
                new = sprintf('e%s%dto%d_', Fnct{i+1}{1}{1}, Fnct{i+1}{2}(1)*1000, Fnct{i+1}{2}(2)*1000);
            else
                new = sprintf('e%dto%d_', Fnct{i+1}{2}(1)*1000, Fnct{i+1}{2}(2)*1000);
            end
        case 'eega_refavg'
            new = 'r_';
        case 'pop_reref'
            new = 'r_';
        case 'eega_normtrltrlvar'
            new = 'nvt_';
        case 'eega_normalization'
            ntyp = 'n';
            ntime = '';
            if any(strcmp(Fnct{i+1},'latency')) && isnumeric(Fnct{i+1}{find(strcmp(Fnct{i+1},'latency'))+1})
                tw = Fnct{i+1}{find(strcmp(Fnct{i+1},'latency'))+1};
                if ismatrix(tw)
                    ntime = sprintf('%dto%d', tw(1), tw(2));
                end
            end
            if any(strcmp(Fnct{i+1},'epochs')) && strcmp(Fnct{i+1}{find(strcmp(Fnct{i+1},'epochs'))+1},'single')
                ntyp = [ntyp 'xep'];
            end
            if any(strcmp(Fnct{i+1},'electrodes')) && strcmp(Fnct{i+1}{find(strcmp(Fnct{i+1},'electrodes'))+1},'single')
                ntyp = [ntyp 'xel'];
            end
            new = [ntyp ntime '_'];
        case 'eega_rmbaseline'
            if any(strcmp('timereference',Fnct{i+1}))
                j = find(strcmp('timereference',Fnct{i+1}));
                new = sprintf('b%s%dto%d_', Fnct{i+1}{j+1}, Fnct{i+1}{1}(1), Fnct{i+1}{1}(2));
            elseif ismatrix(Fnct{i+1}{1})
                new = sprintf('b%dto%d_', Fnct{i+1}{1}(1), Fnct{i+1}{1}(2));
            end
%             if iscellstr(Fnct{i+1}{1})
%                 new = sprintf('b%sto%s_', Fnct{i+1}{1}{1}, Fnct{i+1}{1}{2});
%             elseif ismatrix(Fnct{i+1}{1})
%                 new = sprintf('b%dto%d_', Fnct{i+1}{1}(1), Fnct{i+1}{1}(2));
%             end
        case 'baseline_itt'
            if ismatrix(Fnct{i+1}{1})
                if length(Fnct{i+1})==2 || ~Fnct{i+1}{3}
                    new = sprintf('bj%dto%d_', Fnct{i+1}{1}(1), Fnct{i+1}{1}(2));
                else
                    new = sprintf('bjxall%dto%d_', Fnct{i+1}{1}(1), Fnct{i+1}{1}(2));
                end
            end
        case 'rmbaselinerelative2pre'
            if ismatrix(Fnct{i+1}{1})
                new = sprintf('b%dto%d_', Fnct{i+1}{1}(1), Fnct{i+1}{1}(2));
            end
        case 'dss_denoise_EEG'
            if Fnct{i+1}{2}>=1
                new = sprintf('%d-%d-dss_',Fnct{i+1}{2},Fnct{i+1}{1});
            else
                new = sprintf('%d-%3.2f-dss_',Fnct{i+1}{2},Fnct{i+1}{1});
            end
        case 'eega_spatialsmoothing'
            new = 'ss_';
        case 'eega_runmara'
            new = 'mara_';
        case 'eega_rmcompica'
            new = 'rc_';
        case 'rmvICA'
            new = 'rc_';
        
    end
    fnctname = [new fnctname];
end



end
