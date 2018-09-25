function [res,W,J,X] = instsim(x,model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% initializations and precomputations
tpts = model.vardata.t;
nt = length(tpts);
simsens = model.options.simsens;
X = cell(nt,1);
if simsens
    dX = X;
end
nu = model.vardata.nu;
nc = model.vardata.npool;
np = nu+nc;
nsc = model.data.nh;
u = x(1:nu);
c = x(nu+1:nu+nc);
sc = x(nu+nc+1:end);
if simsens
    [X0,Y0,tms,A,B,C,cY,h,dX0,dY0,sA,sB] = initialmats(u,c,model);
    sX = X;
else
    [X0,Y0,tms,A,B,C,cY,h] = initialmats(u,c,model);
end
nsize = length(X0);
nemu = model.mdvsim.nemu;
nxp = length(model.data);
% default predicted MDVs set to NaN

for i = 1:nt
    X{i} = X0;
    if simsens
        dX{i} = dX0;
    end
    for j = 1:nsize
        X{i}{j}(:,:) = NaN;
        if simsens
            dX{i}{j}(:,:) = NaN;
        end
    end
end





t = 0;
tf = tpts(1);
iter = 0;
fail = false;
i = 1;

while i<=nt && ~fail
    %t = t0;
    if i==1
        tlast=0;
    else
        tlast=tpts(i-1);
    end
    tf = tpts(i);
    if i<nt
        tnext = tpts(i+1);
    else
        tnext = tf;
    end
    iter = 0;
    while t < tf && ~fail
        if simsens
        [t,X0,Y0,h,tms,dX0,dY0] = integstep(t,tf,X0,Y0,h,tms,A,B,C,cY,u,c,tnext,dX0,dY0,sA,sB);
        else
            [t,X0,Y0,h,tms] = integstep(t,tf,X0,Y0,h,tms,A,B,C,cY,u,c,tnext);
        end
        iter = iter+1;
        if h<0.01*(tf-tlast) && iter > 50
            fail = true;
        end
    end
    X{i} = X0;
    if simsens
        dX{i} = dX0;
    end
    i = i+1;
end


% organizing simulated MDVs and sensitivities
for tx = 1:nt
    for i = 1:nsize
        nX = nemu(i);
        nmdv = size(X{tx}{i},2);
        X{tx}{i} = mat2cell(X{tx}{i},ones(1,nX),nmdv);
        X{tx}{i} = X{tx}{i}';
        %X{i} = cell2mat(X{i});
        if simsens
            dX{tx}{i} = dX{tx}{i}';
            dX{tx}{i} = reshape(dX{tx}{i},np,nX*nmdv);
            dX{tx}{i} = mat2cell(dX{tx}{i},np,nmdv*ones(1,nX));
            dX{tx}{i} = [dX{tx}{i}(:)]';
            %sX{i} = cell2mat(sX{i});
        end
    end
end

% collecting data
res = zeros(0,1);
std = res;
if simsens
    J = zeros(0,length(x));
end
ctr = 0;
for xpt = 1:nxp
    nvmeas = length(model.data(xpt).flxind);
    vpred = model.vardata.N*u;
    vpred = vpred(model.data(xpt).flxind);
    vmeas = model.data(xpt).flxval;
    if simsens
        dvpred = [model.vardata.N(model.data(xpt).flxind,:),zeros(nvmeas,sum(nc+nsc))];
    end
    verr = model.data(xpt).flxwt;
    xmeas = model.data(xpt).msval';
    xpred = cell(size(xmeas));
    dxpred = xpred;
    mserr = model.data(xpt).mswt';
    for i = 1:length(xmeas)
        ctr = ctr+1;
        nonan = model.data(xpt).noNaN{i};
        msind = model.data(xpt).msind{i};
        xpred{i} = X{msind(1)}{msind(3),xpt}{msind(4)};
        xpred{i} = conv(xpred{i},model.data(xpt).mscorr{i});
        xk = xpred{i}(1:msind(5));
        xpred{i} = sc(ctr)*xpred{i}(1:msind(5));
        xpred{i} = xpred{i}(nonan);
        mserr{i} = mserr{i}(nonan);
        xmeas{i} = xmeas{i}(nonan);
        if simsens
            dxpred{i} = dX{msind(1)}{msind(3),xpt}{msind(4)};
            dxpred{i} = conv2(dxpred{i},model.data(xpt).mscorr{i});
            dxpred{i} = dxpred{i}(:,1:msind(5));
            dh = zeros(nsc,msind(5));
            dh(ctr,:) = xk;
            dxpred{i} = [sc(ctr)*dxpred{i};dh];
            dxpred{i} = dxpred{i}(:,nonan);
        end
        
    end


    res = [res',(vpred-vmeas)',[cell2mat(xpred)-cell2mat(xmeas)]]';
    std = [std',verr',cell2mat(mserr)]';

    if simsens
         J = [J',dvpred',cell2mat(dxpred)]';
    
    else
        J = 1;
    end
end
W = diag(1./(std.^2));

end

function [X0,Y0,tms,A,B,C,cY,h,dX0,dY0,sF,sG] = initialmats(u,c,model)

simsens = nargout > 7;
A0 = model.mdvsim.A;
B0 = model.mdvsim.B;
A = A0;
B = B0;
sF = cell(size(A));
sG = cell(size(A));
C0 = model.mdvsim.dcxdc;
C = C0;
nb = length(A0);
nu = length(u);
nc = length(c);
cY = model.mdvsim.cY;
X0 = model.mdvsim.X0;
Y0 = model.mdvsim.Y0;
dX0 = model.mdvsim.sX;
dY0 = model.mdvsim.sY;
hinv = 0;   %initial guess set to infinite time step

for i = 1:nb
    nx = size(A0{i},1);
    ny = size(B0{i},2);
    nmdv = length(X0{i}(1,:));
    A0{i} = model.mdvsim.A{i} + sparse(reshape(model.mdvsim.dAdu{i}*u,nx,nx));
    B0{i} = model.mdvsim.B{i} + sparse(reshape(model.mdvsim.dBdu{i}*u,nx,ny));
    [j,k] = find(C0{i}');
    C{i} = sparse(diag(1./c(j)));
    A{i} = C{i}*A0{i};
    B{i} = C{i}*B0{i};
    if simsens
        dX0{i} = zeros(nx,nmdv*(nu+nc));
        dY0{i} = zeros(ny,nmdv*(nu+nc));
        sA = model.mdvsim.sA{i};
        sB = model.mdvsim.sB{i};
        sA = mat2cell(sA,nx*ones(nu,1),nx);
        sA = sA';
        sA = cell2mat(sA);
        sA = C{i}*sA;
        sA = mat2cell(sA,nx,nx*ones(nu,1));
        sA = sA';
        sA = cell2mat(sA);
        sF{i} = sA;
        sB = mat2cell(sB,nx*ones(nu,1),ny);
        sB = sB';
        sB = cell2mat(sB);
        sB = C{i}*sB;
        sB = mat2cell(sB,nx,ny*ones(nu,1));
        sB = sB';
        sB = cell2mat(sB);
        sG{i} = sB;
        sAp = zeros(nx*nc,nx);
        sBp = zeros(nx*nc,ny);
        [j,k] = find(C0{i});
        jn = j + (k-1)*nx;
        sAov = -C{i}*A{i};
        sBov = -C{i}*B{i};
        sAp(jn,:) = sAov(j,:);
        sBp(jn,:) = sBov(j,:);
        sF{i} = [sF{i};sAp];
        sG{i} = [sG{i};sBp];
        %{
        
        for cx = 1:nc
            sfc = -(C{i}*diag(C0{i}(:,cx)))*A{i};
            sF{i} = [sF{i};sfc];
            sgc = -(C{i}*diag(C0{i}(:,cx)))*B{i};
            sG{i} = [sG{i};sgc];
        end
        %}
        Ax = full(A{i} + A{i}')/2;
        hinv = max(hinv,max(abs(eig(Ax))));
        %dX0{i} = [dX0{i},sparse(nx,nc)];
        %dY0{i} = [dY0{i},sparse(ny,nc)];
    end
end

%computing initial time step length
h = 2/hinv;
h = model.vardata.t(1)/ceil(model.vardata.t(1)/h);
%initial transition matrices

tms = tmscalc(A,h,simsens,sF);
    
end

function tms = tmscalc(A,h,simsens,sA)
%warning('off','all');
nsize = length(A);

[tms.F,tms.G,tms.W] = deal(cell(nsize,2));
tms1 = tms;
if simsens
    [tms.sF,tms.sG,tms.sW] = deal(cell(nsize,2));
    np = size(sA{1},1)/size(A{1},1);
end

for i = 1:nsize
    nx = size(A{i},1);
    I = speye(nx);
    Ainv = A{i}\I;
        
    for j = 1:2
    
    % matrix GAMMA and OMEGA
        Exz = expon(A{i}*h/lh);
        tms.F{i,j} = Exz;
        tms.G{i,j} = (tms.F{i,j}-I)*Ainv;
        tms.W{i,j} = ((tms.G{i,j}*j/h)-I)*Ainv;

    end
end
end

function E = expon(A)
% Matrix exponential calculation using Pade's approximation. This code was
% taken directly from expdemo1.m with resquaring from van Loan (1978)
[~,s] = log2(norm(A,inf));
s = max(1,s+1);
invs = 1/(2^s);    %matrix scaling factor to apply Pade's approximation
A = invs*A;

[nx,nx] = size(A);

% Pade approximation for exp(A)


    X = A;

c = 1/2;
E = c*X;
D = -c*X;

E(1:nx,1:nx) = speye(nx) + E(1:nx,1:nx);
D(1:nx,1:nx) = speye(nx) + D(1:nx,1:nx);

q = 6;
p = 1;
for k = 2:q
   c = c * (q-k+1) / (k*(2*q-k+1));
   X = A*X;

   cX = c*X;
   E = E + cX;
   if p
      D = D + cX;
   else
      D = D - cX;
   end
   p = ~p;
end


    E = D\E;
    % Undo scaling by repeated squaring
    for k = 1:s
       E = E*E;
    end
    E(abs(E)<eps) = 0;
    
end


function [tf,X0,Y0,h,tms,dX0,dY0] = integstep(t0,tmax,X0,Y0,h,tms,A,B,C,cY,u,c,tnext,dX0,dY0,sA,sB)

%nu = model.vardata.nu;
%nh = sum([model.data.nh]);
%nc = model.vardata.nc;
%np = nu+nc;
h_old=h;
simsens = nargout>5;
if simsens
    np = length(dX0{1}(1,:))/2;
end
nsize = length(A);

Y1 = Y0;
Y2a = Y0;
Y2b = Y0;
X1s = X0;
X2a = X0;
X2s = X0;
if simsens
    dY1 = dY0;
    dY2a = dY0;
    dY2b = dY0;
    dX1s = dX0;
    dX2a = dX0;
    dX2s = dX0;
end
err = 0;
mdvtol = 1e-3;
stol = 1e-2;
done = false;

% main loop

while ~ done
    Ediff1 = zeros(0,1);
    if simsens
        Ediff1 = zeros(0,1+np);
    end
    
    for i = 1:nsize
        
        %protect against negative mass fractions
        X0{i} = max(X0{i},0);
        
        nmdv = length(X0{i}(1,:));
        nX = length(X0{i}(:,1));
        %single step
        F1 = tms.F{i,1};
        G1 = tms.G{i,1};
        W1 = tms.W{i,1};
        F2 = tms.F{i,2};
        G2 = tms.G{i,2};
        W2 = tms.W{i,2};
        %{
        if simsens
            dF1 = tms.sF{i,1};
            dG1 = tms.sG{i,1};
            dW1 = tms.sW{i,1};
            dF2 = tms.sF{i,2};
            dG2 = tms.sG{i,2};
            dW2 = tms.sW{i,2};
        end
%}
        X1s{i} = F1*X0{i} + G1*B{i}*Y0{i} + W1*B{i}*(Y1{i}-Y0{i});
        % Ensure that MDVs always sum to unity
        %X1s{i} = max(X1s{i},0);
        Sm = sum(X1s{i},2);
        Sx = diag(Sm)*ones(length(Sm),length(X1s{i}(1,:)));
        X1s{i} = X1s{i}./Sx;
        if simsens
            %{
            D = dF1*X0{i} + dG1*B{i}*Y0{i} + dW1*B{i}*(Y1{i}-Y0{i});
            D = reshape(D,nX,nmdv*np);
            D1 = sB{i}*Y0{i};
            D2 = sB{i}*(Y1{i}-Y0{i});
            D1 = G1*reshape(D1,nX,nmdv*np);
            D2 = W1*reshape(D2,nX,nmdv*np);
            dX1s{i} = D + D1 + D2 + F1*(dX0{i}) + G1*B{i}*(dY0{i}) + W1*B{i}*(dY1{i}-dY0{i});
            %}
            %
            %H0 = sA{i}*X0{i} + sB{i}*Y0{i};
            %H0 = reshape(H0,nX,nmdv*np);
            %H0 = H0 + B{i}*dY0{i};
            %H1 = sA{i}*X1s{i} + sB{i}*Y1{i};
            %H1 = reshape(H1,nX,nmdv*np);
            %H1 = H1 + B{i}*dY1{i};
            H0 = reshape(sA{i}*X1s{i} + sB{i}*Y1{i},nX,nmdv*np);
            H1 = reshape(sA{i}*X0{i} + sB{i}*Y0{i},nX,nmdv*np);
            H = G1*H1 + W1*(H0-H1);
            dX1s{i} = F1*dX0{i} + G1*B{i}*dY0{i} + W1*B{i}*(dY1{i}-dY0{i}) + H;
            %}
            
        end
        if i<nsize
            if simsens
                [Y1{i+1},dY1{i+1}] = emuconv(X1s,i+1,cY{i+1},Y1{i+1},dX1s,dY1{i+1});
            else
                Y1{i+1} = emuconv(X1s,i+1,cY{i+1},Y1{i+1});
            end
        end

        %two steps

        %first step
        X2a{i} = F2*X0{i} + G2*B{i}*Y0{i} + W2*B{i}*(Y2a{i}-Y0{i});
        %X2a{i} = max(X2a{i},0);
        Sm = sum(X2a{i},2);
        Sx = diag(Sm)*ones(length(Sm),length(X2a{i}(1,:)));
        X2a{i} = X2a{i}./Sx;
        if simsens

            H0 = reshape(sA{i}*X2a{i} + sB{i}*Y2a{i},nX,nmdv*np);
            H1 = reshape(sA{i}*X0{i} + sB{i}*Y0{i},nX,nmdv*np);
            H = G2*H1 + W2*(H0-H1);
            dX2a{i} = F2*dX0{i} + G2*B{i}*dY0{i} + W2*B{i}*(dY2a{i}-dY0{i}) + H;
            %}
        end
        if i<nsize
            if simsens
                [Y2a{i+1},dY2a{i+1}] = emuconv(X2a,i+1,cY{i+1},Y2a{i+1},dX2a,dY2a{i+1});
            else
                Y2a{i+1} = emuconv(X2a,i+1,cY{i+1},Y2a{i+1});
            end
        end

        %second step
        X2a{i} = max(X2a{i},0);
        X2s{i} = F2*X2a{i} + G2*B{i}*Y2a{i} + W2*B{i}*(Y2b{i}-Y2a{i});
        %X2s{i} = max(X2s{i},0);
        Sm = sum(X2s{i},2);
        Sx = diag(Sm)*ones(length(Sm),length(X2s{i}(1,:)));
        X2s{i} = X2s{i}./Sx;

        if simsens

            H0 = reshape(sA{i}*X2s{i} + sB{i}*Y2b{i},nX,nmdv*np);
            H1 = reshape(sA{i}*X2a{i} + sB{i}*Y2a{i},nX,nmdv*np);
            H = G2*H1 + W2*(H0-H1);
            dX2s{i} = F2*dX2a{i} + G2*B{i}*dY2a{i} + W2*B{i}*(dY2b{i}-dY2a{i}) + H;
            %}
        end
        if i<nsize
            if simsens
                [Y2b{i+1},dY2b{i+1}] = emuconv(X2s,i+1,cY{i+1},Y2b{i+1},dX2s,dY2b{i+1});
            else
                Y2b{i+1} = emuconv(X2s,i+1,cY{i+1},Y2b{i+1});
            end
        end

        %error calculation
        Ediff = (X2s{i}-X1s{i})/mdvtol;
        wt = 1./max(abs(X0{i}),1);
        Ediff = Ediff(:);
        wt = wt(:);
        if simsens
            Sx2s = dX2s{i}-dX1s{i};
            %Sx1s = dX1s;
            Sx0 = dX0{i};
            Sx2s = mat2cell(Sx2s,ones(1,nX),np*ones(1,nmdv));
            Sx2s = cell2mat(Sx2s(:));
            Sx2s = Sx2s*diag([u;c]);
            Sx0 = mat2cell(Sx0,ones(1,nX),np*ones(1,nmdv));
            Sx0 = cell2mat(Sx0(:))*diag([u;c]);
            Sx0 = 1./max(abs(Sx0),1);
            Ediff = [Ediff,Sx2s/stol];
            wt = [wt,Sx0];
        end
        Ediff = Ediff.*wt/3;
        Ediff1 = [Ediff1;Ediff];
        %e1 = rms(Ediff,1);
        %err = max(err,max(e1));
    end
    
    %error calculation
    e1 = rms(Ediff1,1);
    err = max(e1);

    %step rescaling

    if err <= 1
        tf = t0+h;
        if tf < tmax
            hmax = (tmax-tf);
        else
            hmax = tnext-tf;
        end
            sc = min(1/sqrt(err),1000);
            h1 = h*sc;
            ha = min(0.9*h1,hmax);
            if hmax > 0
                ha = hmax/(ceil(hmax/ha));
            else
                ha = 0;
            end
            kx = hmax/h;
            %if mod(hmax/h,1)>0
            if abs(kx-round(kx))>sqrt(eps)*tf
                hb = ha;
            else
                hb = h;
            end
            if h1 >= 2.2*h
                h = ha;
            else
                h = hb;
            end
            if (hmax > 0) && ~isequal(h_old,h)
                if simsens
                    tms = tmscalc(A,h,simsens,sA);
                else
                    tms = tmscalc(A,h,0,0);
                end
            end
                

            

        X0 = X2s;
        Y0 = Y2b;
        if simsens
            dX0 = dX2s;
            dY0 = dY2b;
        end
        done = true;

    else
        sc = min(max(1/sqrt(err),0.001),0.5);
        err=0;
        h = h*sc;
        if simsens
            tms = tmscalc(A,h,simsens,sA);
        else
            tms = tmscalc(A,h,0,0);
        end
        done = false;
    end

end
end

function [y,sy] = emuconv(X,esize,C,y,sX,sy)
% helper function to convolute EMUs of smaller sizes to serve as inputs to
% the EMU network of larger size

if nargin<5
    simsens = false;
else
    simsens = true;
end

%X = X(1:(esize-1));
if simsens
    %sX = sX(1:(esize-1));
    np = length(sX{1}(1,:))/2;
end
for i = 1:(esize-1)
    [nx,nmdv] = size(X{i});
    X{i} = mat2cell(X{i},ones(1,nx),nmdv);
    X{i} = X{i}';
    if simsens
        sX{i} = sX{i}';
        sX{i} = reshape(sX{i},np,nx*nmdv);
        sX{i} = mat2cell(sX{i},np,nmdv*ones(1,nx));
        %sX{i} = cell2mat(sX{i}(:)');
        %sX{i} = mat2cell(sX{i},1,(i+1)*ones(1,nx));
    end
end
%X = cell2mat(X');
X1 = X;
if simsens
sX1 = sX;
end
X = cell(1,0);
sX = cell(1,0);
for i = 1:esize-1
    X = [X,X1{i}];
    if simsens
        sX = [sX,sX1{i}];
    end
    %X{i} = [];
end
nmdv = length(X1{esize}(1,:));
%convoluting EMUs
ny = size(C,1);
%y = zeros(ny,length(X1{esize}(1,:)));
%if simsens
%    sy = zeros(ny,np*(esize+1));
%end
inds = 1:ny;
inds = inds';
i1 = inds(any(C,2));
for i = 1:length(i1)
    ix = i1(i);
    in2 = find(C(ix,:));
    ncv = C(ix,in2);
    i2 = zeros(1,sum(ncv));
    cv = 0;
    for j = 1:length(in2)
        i2(cv+1:cv+ncv(j)) = in2(j);
        cv = cv+ncv(j);
    end
    y0 = 1;
    if simsens
        sy0 = zeros(np,nmdv);
    end
    for j = 1:length(i2)
        y0 = conv(y0,X{i2(j)});
        if simsens
            syk = 1;
            for k = 1:length(i2)
                if k == j
                    syk = conv2(syk,sX{i2(k)});
                else
                    syk = conv2(syk,X{i2(k)});
                end
            end
            sy0 = sy0 + syk;
        end
    end
    y(ix,:) = y0;
    if simsens
        sy(ix,:) = sy0(:)';
    end
end
%y = cell2mat(y);
end


