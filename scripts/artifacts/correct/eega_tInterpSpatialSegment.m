% -------------------------------------------------------------------------
% The function interpolates for each channels the bad segments if they are
% comprise in the InterTime and if there are enought good channels for
% interpolation
%
% INPUT
% DATA      :
% BCT : bad data
% GoodCh :
% InterTime :
% chanlocs :
% p_int : maximun proportion of bad channels to interpolate
%
% OPTIONAL INPUTS
% distmethod    : methos used to compute the neighbours ('triangulation' or
%                 'distance', default 'triangulation')
% pneigh        : maximun proportion of bad neighbor channels
% distneighbour : distance to define neighbours
% minintertime  : minimun length of a segment in order to be interpolated
%               in samples (default=1).
% masktime      : time to mask bad segments before interpolation in samples
%                 (default=0)
% maxloop       : maximun number of loops for interpolation(default=10)
% splicemethod  : how the interpolated segments segments are splice in the
%                data
%               - 0 = none
%               - 1 = segments are aligened with the previous segment
%               - 2 = segments are aligned in the midle of the previous and the next
%               - 3 = segments are linearly fit between the two
%
% OUTPUT:
% DATA      : data matrix with the interpolated data
% Interp    : logical matrix indicating the data that was interpolated
%
% USAGE
% [ DATA, Interp ] = eega_tInterpSpatialSegment(DATA, BCT, GoodCh, InterTime, chanlocs, 0.4, 'distmethod', 'triangulation')
%
% -------------------------------------------------------------------------

function [ DATAgood, Interp ] = eega_tInterpSpatialSegment(DATA, BCT, GoodCh, InterTime, chanlocs, p_int, varargin )

fprintf('### Spatial interpolatiom of bad segments ###\n')

%% ------------------------------------------------------------------------
%% Parameters

P.minintertime      = 0;
P.masktime          = 0;
P.distneighbour     = [];
P.distmethod        = 'triangulation'; % 'distance', 'triangulation' or 'template'
P.pneigh            = 1;
P.maxloop           = 10;
P.splicemethod      = 1; % 0 / 1 / 2 / 3
P.silent            = 0;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tInterpSpatialSegment: Non recognized inputs')
end

% Check the input
if ~any(P.splicemethod==[0 1 2 3 4])
    error('eega_tInterpSpatialSegment: wrong splicing method. Options 0 | 1 | 2| 3 | 4')
end
if strcmp(P.distmethod,{'distance'  'triangulation'  'template' })
    error('eega_tInterpSpatialSegment: distmethod can have values ''distance'' / ''triangulation'' / ''template'' ')
end

%% ------------------------------------------------------------------------
%% Obtain electrodes position and G
Xel = [chanlocs.X];
Yel = [chanlocs.Y];
Zel = [chanlocs.Z];
rad = sqrt(Xel.^2+Yel.^2+Zel.^2);
Xel =  Xel./rad;
Yel =  Yel./rad;
Zel =  Zel./rad;
Gel = computeg(Xel,Yel,Zel,Xel,Yel,Zel);

%% ------------------------------------------------------------------------
%% Define neighbours

nEl = size(DATA,1);
nS = size(DATA,2);
nEp = size(DATA,3);

if ~isempty(P.pneigh) && (P.pneigh<1)
    ChLabels = cell(nEl,1);
    Chanpos = cat(2,Xel(:),Yel(:),Zel(:));
    for el=1:nEl
        ChLabels{el} = chanlocs(el).labels;
    end
    layout.label = ChLabels;
    layout.chanpos = Chanpos;
    neighbours = eega_neighborchannels(layout,'neighbourdist',P.distneighbour,'method',P.distmethod);
else
    neighbours = [];
end


%% ------------------------------------------------------------------------
%% Interpolation

ok = 0;
Interp = false([nEl nS nEp]);
datainterp = sum(Interp(:));
counter = 0;
DATAgood = DATA;
while ~ok  && counter<P.maxloop %loop till all possible data to interpolate is interpolated
    counter = counter+1;
    
    %% Find the segments to interpolate
    
    bctint = false([nEl nS nEp]);
    if ~P.silent
        fprintf('Total bad time % 8d (% 7.2f%%).\n',sum(sum(~InterTime(1,:,:),2),3),sum(sum(~InterTime(1,:,:),2),3)/(nS*nEp)*100)
    end
    for el=1:nEl
        
        for ep = 1:nEp
            tointer = BCT(el,:,ep) & InterTime(1,:,ep) & repmat(GoodCh(el,:,ep),[1 nS 1]) & ((sum(BCT(:,:,ep),1)/nEl)<=p_int);
            
            %begining end of each segment
            badi = find( [ tointer(1) diff(tointer)==1 ] )';
            badf = find( [ diff(tointer)==-1 tointer(end) ] )';
            
            %remove too short segments
            baddur_ep = badf - badi + 1;
            idrmv = baddur_ep < P.minintertime;
            badi(idrmv,:) = [];
            badf(idrmv,:) = [];
            
            %mask
            badif_ep = [badi-P.masktime badf+P.masktime];
            badif_ep(badif_ep(:,1)<1,1) = 1;
            badif_ep(badif_ep(:,2)>nS,2) = nS;
            
            if ~isempty(badif_ep)
                
                % put together segments which overlap
                if size(badif_ep,1)>1
                    id = find(badif_ep(1:end-1,2)>=badif_ep(2:end,1));
                    for iii=1:length(id)
                        badif_ep(id(iii),2) = badif_ep(id(iii)+1,2);
                    end
                    badif_ep(id+1,:) = [];
                end
                
                % mark in the matrix for interpolation
                for i=1:size(badif_ep,1)
                    bctint(el,badif_ep(i,1):badif_ep(i,2),ep) = 1 ;
                end
                
            end
            
        end % epoch
        
        if ~P.silent
            fprintf('Electrode % 4.0d. Bad data % 8d (% 7.2f%%). Data to interpolate % 8d (% 7.2f%%)\n',...
                el,...
                sum(sum(BCT(el,:,:),2),3), sum(sum(BCT(el,:,:),2),3)/(nS*nEp)*100,...
                sum(sum(bctint(el,:,:),2),3), sum(sum(bctint(el,:,:),2),3)/(nS*nEp)*100);
        end
        
    end % electrode
    
    %% Interpolate
    % Interp = zeros([nEl nS nEp]);
    for ep = 1:nEp
        t0=tic;
        fprintf('Interpolating epoch %d...',ep)
        
        Iall = nan(1,nS);
        Fall = nan(1,nS);
        thesement = 0;
        
        % Indentify changes in the number of bad channels
        % -----------------------------------------------
        change_ch = find(any(diff(bctint(:,:,ep),1,2),1));
        seg_i = unique([change_ch+1 1],'sorted');
        seg_f = unique([change_ch nS],'sorted');
        bad_if = [seg_i' seg_f'];
        
        % interpolate
        % ----------------------
        for iseg=1:size(bad_if,1)
            if any(any(bctint(:,bad_if(iseg,1):bad_if(iseg,2),ep),2),1)
                
                thesement = thesement+1;
                
                Im = bad_if(iseg,1);
                Fm = bad_if(iseg,2);
                Iall(thesement) = Im;
                Fall(thesement) = Fm;
                
                AllBadCh    = any(bctint(:,Im:Fm,ep),2) | ~GoodCh(:,1,ep);
                BadCh2Int   = any(bctint(:,Im:Fm,ep),2) & GoodCh(:,1,ep);
                
                [badchandata, BadCh2Int, interpdone] = eega_InterpSphericSpline(DATA(:,Im:Fm,ep), BadCh2Int, AllBadCh, Xel, Yel, Zel, 'G', Gel, 'p_int', p_int, 'pneigh', P.pneigh, 'neighbours',neighbours);
                
                if interpdone
                    Interp(BadCh2Int,Im:Fm,ep) = 1;
                    DATAgood(BadCh2Int,Im:Fm,ep) = badchandata;
                    BCT(BadCh2Int,Im:Fm,ep) = 0;
                end
            end
        end % bad segments
        
        Iall(isnan(Iall)) = [];
        Fall(isnan(Fall)) = [];
        
        % splice the interpolated segments together
        % -----------------------------------------
        switch P.splicemethod
            case 1 % align to the previous
                DATAgood(:,:,ep) = eega_splicesgments1(DATAgood(:,:,ep), [Iall(:) Fall(:)], [1 size(DATA,2)], BCT(:,:,ep), ~InterTime(:,:,ep), ~GoodCh(:,:,ep));
            case 2 % align in the midle of the previous and the next
                DATAgood(:,:,ep) = eega_splicesgments2(DATAgood(:,:,ep), [Iall(:) Fall(:)], [1 size(DATA,2)], BCT(:,:,ep), ~InterTime(:,:,ep), ~GoodCh(:,:,ep));
            case 3 % linearly fit between the two
                DATAgood(:,:,ep) = eega_splicesgments3(DATAgood(:,:,ep), [Iall(:) Fall(:)], [1 size(DATA,2)], BCT(:,:,ep), ~InterTime(:,:,ep), ~GoodCh(:,:,ep));
        end
        
        fprintf('\n')
        t=toc(t0);
        fprintf('--> interpolation time %f\n',t)
        fprintf('--> percentage of interpolated data %3.1f %%\n',sum(sum(Interp(:,:,ep)))/(nEl*nS)*100)
        fprintf('\n')
    end % epochs
    
    %% Check waht was interpolated
    BCT(Interp & BCT) = 0;
    datainterpnew = sum(Interp(:));
    if ((datainterpnew - datainterp) / numel(Interp) *100) == 0
        ok = 1;
    else
        datainterp = datainterpnew;
    end
    
end
end

