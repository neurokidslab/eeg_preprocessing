% -------------------------------------------------------------------------
% This function splices segments of data at the indicated point by aligning
% each segemtn by lienarly twisting them to be continuisly concatenated
% with the segment before and the segment after
%
% INPUT
% data      : data (electrodes x samples)
% bad_if    : matrix indicating the begingn and end of the segments 
%             (n segments x 2)
% epoch_if  : matrix indicating the begingn and end of the epochs
% bct       : matrix indicating bad data (electrodes x samples)
%             If it is empty is not used
% bt        : matrix indicating bad time (1 x samples)
%             If it is empty is not used
% bc        : matrix indicating bad electrodes (electrodes x 1)
%             If it is empty is not used
%
% OUTPUT
% dN        : data corrected
%
% -------------------------------------------------------------------------

function dN = eega_splicesgments3(data, bad_if, epoch_if, bct, bt, bc)
if nargin<4 || isempty(bct)
    bct = false(size(data));
end
if nargin<5 || isempty(bt)
    bt = false(1, size(data, 2));
end
if nargin<6 || isempty(bc)
    bc = false(size(data, 1), 1);
end

Iall = bad_if(:,1)';
Fall = bad_if(:,2)';

epochI = epoch_if(:,1)';
epochF = epoch_if(:,2)';

% splice segments
dN = data;
if ~isempty(Iall)
    for i=1:length(Iall)
        yadd=nan(size(data,1), 2);
        % take the values before
        if ~any(ismember(epochI,Iall(i)))
            ide = ~bct(:,Iall(i)-1) & ~bt(Iall(i)-1) & ~bc; % only for the electrodes with good data
            if ~isempty(ide)
                yadd(ide, 1) = data(ide, Iall(i)-1);
            end
        end
        % take the values after
        if ~any(ismember(epochF,Fall(i)))
            ide = ~bct(:,Fall(i)+1) & ~bt(Fall(i)+1) & ~bc; % only for the electrodes with good data
            if ~isempty(ide)
                yadd(ide, 2) = data(ide, Fall(i)+1);
            end
        end
        % align the beging of the bad segment with the end of the previous segment
        ide = ~isnan(yadd(:,1));
        if ~isempty(ide)
            dN(ide, Iall(i):Fall(i)) = data(ide, Iall(i):Fall(i)) - data(ide, Iall(i)) + yadd(ide, 1);
        end
        % align the end to the beging of the next segment
        ide = ~isnan(yadd(:,2));
        if ~isempty(ide)
            p =  ( yadd(ide, 2) - dN(ide, Fall(i)) ) / (Fall(i) - Iall(i) + 1);
            dN(ide,Iall(i):Fall(i)) = dN(ide, Iall(i):Fall(i)) + p .* (1:(Fall(i) - Iall(i) + 1));
        end
        
    end
end

end