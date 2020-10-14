function [P, OK, extrainput] = eega_getoptions(P, inputs)

if mod(length(inputs),2)==1
    error('eega_getoptions: Optional parameters come by pairs')
end

fP = fieldnames(P);
fV = inputs(1:2:end);
vV = inputs(2:2:end);
fVfP = any(strcmpi(repmat(fV(:),[1 length(fP)]),repmat(fP(:)',[length(fV(:)) 1])),2);
Pop = [];
for i=1:length(fV)
    Pop.(fV{i}) = vV{i};
end
extrainput = {};
j=1;
if ~all(fVfP)
    badinput = find(~fVfP);
    for i=1:length(badinput)
        extrainput{j} = inputs{badinput(i)*2-1};
        extrainput{j+1} = inputs{badinput(i)*2};
        j=j+2;
    end
    OK = 0; 
else
    OK = 1;
end
for f=1:numel(fP)
    if isfield(Pop,fP{f})
        P.(fP{f}) = Pop.(fP{f});
    end
end

end