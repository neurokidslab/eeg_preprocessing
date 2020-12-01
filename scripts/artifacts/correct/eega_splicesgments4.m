% -------------------------------------------------------------------------
% This function splices segments of data at the indicated point by first 
% roboust detrending the data using good and not interpolated segments, and
% then linearly fitting the interpolated segment within the data 
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
% nplus     : number of samples to mask the bad segment to detrend
% order     : orde of polynomial to fit the trend
%
% OUTPUT
% dN        : data corrected
%
% -------------------------------------------------------------------------

function dN = eega_splicesgments4(data, bad_if, epoch_if, bct, bt, bc, window, order)

% check that the NoiseTools are in the path
if exist('nt_detrend','file')~=2
    error('The NoiseTools is not in the path. Add them to the path (http://audition.ens.fr/adc/NoiseTools/)')
end

% parameters
if nargin<4 || isempty(bct)
    bct = false(size(data));
end
if nargin<5 || isempty(bt)
    bt = false(1, size(data, 2));
end
if nargin<6 || isempty(bc)
    bc = false(size(data, 1), 1);
end
if nargin<7 || isempty(window)
    window = 4000;
end
if nargin<8 || isempty(order)
    order = 3;
end
window = floor(window/2)*2;



Iall = bad_if(:,1)';
Fall = bad_if(:,2)';

epochI = epoch_if(:,1)';
epochF = epoch_if(:,2)';

dN = data;
szd = size(data);
if ~isempty(Iall)
    for i=1:length(Iall)
        
        % TAKE A SOURRONDING SEGMENT
        % --------------------------
        nbad = (Fall(i) - Iall(i)) +1;
        nplus = window/2;
        w = [];
        x = [];
        % take the values before
        if ~any(ismember(epochI,Iall(i)))
            s1 = max(Iall(i)-nplus,1);
            w = cat(2,w,false(szd(1), Iall(i)-s1));
            x = cat(2,x,data(:,s1:Iall(i)-1));
        else
            s1 = Iall(i);
        end
        % take the values during the bad segment
        w = cat(2,w,false(szd(1),nbad));
        x = cat(2,x,data(:,Iall(i):Fall(i)));
        % take the values after
        if ~any(ismember(epochF,Fall(i)))
            s2 = min(Fall(i)+nplus,szd(2));
            w = cat(2,w,false(szd(1), s2-Fall(i)));
            x = cat(2,x,data(:,Fall(i)+1:s2));
        else
            s2 = Fall(i);
        end
        bd = bct(:,s1:s2) | repmat(bt(s1:s2),[szd(1) 1]) | repmat(bc, [1 (s2-s1+1)]);
        sx1 = Iall(i)-s1+1;
        sx2 = size(x,2) - (s2-Fall(i)) ;
        
        % SPLICE THE BAD SEGMENT BY CONCATENATING
        % ---------------------------------------
        y = eega_splicesgments1(x, [sx1 sx2], [1 size(x,2)], bct(:,s1:s2), bt(:,s1:s2), bc);
        
        % ROBOUST DETRENDING OF THE BIG SEGMENT
        % -------------------------------------
        [z,wout,~] = nt_detrend(y', order, ~bd');
        z = z';
        
        % LINEARLY FIT THE SEGMENT IN THE REST OF THE DATA
        % ------------------------------------------------
        dN(:,s1:s2) = z;
        dN = eega_splicesgments3(dN, [s1 s2], epoch_if, bct, bt, bc);
        
    end
end


% % detrend good parts
% w = true(size(data));
% for i=1:length(Iall)
%     w(:,Iall(i):Fall(i)) = 0;
% end
% w(bct) = 0;
% w(:,bt) = 0;
% w(bc,:) = 0;
% [y,wout] = nt_detrend(data',order,w',[],[],[],window);
% y = y';
% wout = wout';
% y(wout==0) = data(wout==0);
% 
% % find the segments that were not detrended
% sy1 = find(diff([0 ~wout])==1);
% sy2 = find(diff([~wout 0])==-1);
%          
% % align the bad segments within the sourronding ones
% z = eega_splicesgments3(y, [sy1(:) sy2(:)], epoch_if, bct, bt, bc);
% dN = z;


% % splice segments
% dN = data;
% szd = size(data);
% if ~isempty(Iall)
%     for i=1:length(Iall)
%         
%         % DETREND THE PERIOD AROUND THE BAD SEGMENT
%         nbad = (Fall(i) - Iall(i)) +1;
% %         nplus = nbad; % take for detrending a a segment which is 3 times the bad segment
%         w = [];
%         x = [];
%         % take the values before
%         if ~any(ismember(epochI,Iall(i)))
%             s1 = max(Iall(i)-nplus,1);
%             nw = ~bct(:,s1:Iall(i)-1) & ~bt(s1:Iall(i)-1) & ~bc;
%             w = cat(2,w,nw);
%             x = cat(2,x,data(:,s1:Iall(i)-1));
%         else
%             s1 = Iall(i);
%         end
%         % take the values during the bad segment
%         nw = false(szd(1),nbad);
%         w = cat(2,w,nw);
%         x = cat(2,x,data(:,Iall(i):Fall(i)));
%         % take the values after
%         if ~any(ismember(epochF,Fall(i)))
%             s2 = min(Fall(i)+nplus,szd(2));
%             nw = ~bct(:,Fall(i)+1:s2) & ~bt(Fall(i)+1:s2) & ~bc;
%             w = cat(2,w,nw);
%             x = cat(2,x,data(:,Fall(i)+1:s2));
%         else
%             s2 = Fall(i);
%         end
%         [y,wout,~] = nt_detrend(x',order,w');
%         y = y'; 
%         wout = wout';
%         y(wout==0) = x(wout==0);
%         sy1 = find(diff([0 ~wout])==1);
%         sy2 = find(diff([~wout 0])==-1);
%            
%         % ALIGN THE BAD SEGMENT WITHIN THE BIGGER ONE
%         z = eega_splicesgments3(y, [sy1(:) sy2(:)], [1 size(y,2)], bct(:,s1:s2), bt(:,s1:s2), bc);
%         dN(:,s1:s2) = z;
%         
%         % ALIGN THE DETRENDED SEGMENT WITHIN THE WHOLE DATA
%         dN = eega_splicesgments3(dN, [s1 s2], epoch_if, bct, bt, bc);
%                
%     end
% end

end