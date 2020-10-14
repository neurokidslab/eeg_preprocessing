% -------------------------------------------------------------------------
% This function splices segments of data at the indicated point by aligning
% with the previous segment
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

function dN = eega_splicesgments1(data, bad_if, epoch_if, bct, bt, bc)
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
    for ep=1:length(epochI)
        smplI = Iall((Iall>=epochI(ep)) & (Iall<=epochF(ep)));
        smplF = Fall((Fall>=epochI(ep)) & (Fall<=epochF(ep)));
        theI = unique( [ epochI(ep) smplI smplF+1 epochF(ep)+1] ) ;
        
        if length(theI)>2
            for kk=2:length(theI)-1
                % identify which electrodes are fine and which ones are bad
                if bt(theI(kk)-1)     % if the segment before is bad all electrodes are bad
                    idEl_g = false(size(dN,1),1);
                    idEl_b = true(size(dN,1),1);
                else
                    idEl_g = ~bct(:,theI(kk)-1);    % electrodes with good values before
                    idEl_b = bct(:,theI(kk)-1);     % electrodes with bad values before
                end
                % alignt the good electrodes with the previous
                if any(idEl_g & ~bc)  
                    el = idEl_g & ~bc;
                    yadd = dN(el, theI(kk)-1) - data(el, theI(kk));
                    dN(el,theI(kk):theI(kk+1)-1) = data(el,theI(kk):theI(kk+1)-1) + repmat(yadd,[1 theI(kk+1)-theI(kk)]);
                end
                % if it is a bad electorde align to the good segment before and after
                if any(idEl_b & ~bc)  
                    el = find(idEl_b & ~bc);
                    x = theI(kk-1):theI(kk)-1;
                    if x(1)==1 && x(end)~=size(dN,2)
                        yf = dN(el, x(end)+1);
                        dN(el, x) = dN(el,x) + repmat(yf - dN(el, x(end)), [1 x(end)]);
                    elseif x(end)==size(dN,2) && x(1)~=1
                        yf = dN(el, x(1)-1);
                        dN(el, x) = dN(el,x) + repmat(yf - dN(el, x(1)), [1 x(end)]);
                    elseif x(end)~=size(dN,2) && x(1)~=1
                        yf = mean([dN(el,x(1)-1), dN(el, x(end)+1)], 2);
                        dN(el, 1:x(end)) = dN(el,1:x(end)) + repmat(yf - dN(el,x(1)-1), [1 x(end)]);
                        dN(el, x(end)+1:end) = dN(el, x(end)+1:end) + repmat(yf - dN(el,x(end)+1), [1 size(dN,2)-x(end)]);
                        p = (dN(el,x(end)+1) - dN(el,x(end))) / (length(x)+1);
                        yadd = repmat(p, [1 length(x)]) .* repmat((x - x(1) + 1), [length(p) 1]); 
                        dN(el,x) = dN(el,x) + yadd;
                    end
                end
            end 
        end
    end
end

end