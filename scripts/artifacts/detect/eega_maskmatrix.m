function [BCTm, mask] = eega_maskmatrix(BCT, tmask, srate)

[nEl, nS, nEp] = size(BCT);
mask = false(nEl,nS,nEp);

art_buffer = round(tmask*srate); % time in seconds times sample rate
for ep=1:nEp
    for el = 1:nEl
        bad_i = diff([0; BCT(el,:,ep)'])==1;  % begining of bad segments
        bad_f = diff([BCT(el,:,ep)'; 0])==-1; % end of bad segments
        bad_i = find(bad_i);
        bad_f = find(bad_f);
        for j=1:length(bad_i)
            bad_idx_i = ((bad_i(j)-art_buffer):(bad_i(j)-1));
            bad_idx_i(bad_idx_i<=0) = [];
            bad_idx_f = ((bad_f(j)+1):(bad_f(j)+art_buffer));
            bad_idx_f(bad_idx_f>nS) = [];
            mask(el,bad_idx_i,ep) = 1;
            mask(el,bad_idx_f,ep) = 1;
        end
    end
end
BCTm = BCT | mask;

end