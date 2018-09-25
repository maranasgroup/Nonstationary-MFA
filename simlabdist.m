function [ r,W,J,X ] = simlabdist( x,model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if model.options.ss
    if model.options.simsens
    [r,W,J,X] = stsim(x,model);
    else
        [r,W,~,X] = stsim(x,model);
        J = 0;
    end
else
    if model.options.simsens
    [r,W,J,X] = instsim(x,model);
    else
        [r,W,~,X] = instsim(x,model);
        J = 0;
    end
end

%[r,J] = measextract(X,dXdp.data);

end