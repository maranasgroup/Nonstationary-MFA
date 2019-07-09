function [res,impres] = confintestimate(res,model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% determining the set of reactions for which confidence intervals must be
% computed

N = model.vardata.N;
vfwd = model.vardata.vfwd;
vrev = model.vardata.vrev;
N1 = N;
N1(vfwd,:) = N(vfwd,:)-N(vrev,:);
m = [model.vardata.flxdata.main]';
dil = [model.vardata.flxdata.dilution]';
defid = false(size(vfwd));
set = model.options.conf_set;
if ismember(set,'all')
    cid = m|dil;
    Aex = sparse(diag(cid));
elseif ismember(set,'main')
    cid = m;
    Aex = sparse(diag(cid));
elseif ismember(set,'dilution')
    cid = dil;
    Aex = sparse(diag(cid));
elseif ismember(set,'all_net')
    cid = m;
    cid(vrev) = false;
    Aex = sparse(diag(cid));
elseif ismember(set,'all_exch')
    cid = vrev;
    Aex = sparse(diag(cid));
elseif ismember(set,'minset_main')
    %vconf = minconfset(model);
    cid = model.minset;
    %cid(vconf) = true;
    cid(:,dil) = 0;
    Aex = sparse(cid);
elseif ismember(set,'minset_all')
    %vconf = minconfset(model);
    cid = model.minset;
    Aex = sparse(cid);
else
    % custom list
    customlist = model.options.conf_custom;
    cid = defid;
    cid(customlist) = true;
    Aex = sparse(diag(cid));
end
a1 = any(Aex,2);
Aex = Aex(a1,:);
Aex(:,vrev) = Aex(:,vrev) - Aex(:,vfwd);

% handling pool size ranges for instationary MFA

% setting up constraint matrices

A = [N;-N];
b = [model.vardata.vb(:,1);-model.vardata.vb(:,2)];
nh = length(res.reinit_data.h);

if ~model.options.ss
    c = res.reinit_data.c;
    nc = length(c);
    A = blkdiag(A,eye(nc),eye(nh));
    b = [b;1e-7*ones(nc+nh,1)];
    N1 = [N1,zeros(length(N1(:,1)),nc+nh)];
else
    c = [];
    A = blkdiag(A,eye(nh));
    b = [b;1e-7*ones(nh,1)];
    N1 = [N1,zeros(length(N1(:,1)),nh)];
end
x = [res.reinit_data.u;c;res.reinit_data.h];
actcon = find(A*x<=b);

% estimating preliminary bounds for confidence intervals

%confs = find(cid);
nconfs = length(Aex(:,1));
S_bal=model.vardata.S_bal;
S_balred=S_bal(:,~vrev);
gmod.A = sparse(S_balred);
gmod.rhs = zeros(size(gmod.A(:,1)));
gmod.sense(1:length(gmod.rhs)) = '=';
gmod.lb = [res.fluxes(~vrev).vLB]';
gmod.ub = [res.fluxes(~vrev).vUB]';
gmod.obj = zeros(size(gmod.lb));
%Aex = zeros(nconfs,length(gmod.lb));
gmod.vtype(1:length(gmod.lb)) = 'C';
params.outputflag = 0;
yb = zeros(nconfs,2);
for i = 1:nconfs
    g1 = gmod;
    g1.obj = full(Aex(i,~vrev))';
    g1.obj=double(g1.obj);
    %if vfwd(confs(i)) == true
    %    g1.obj(confs(i)+1) = -1;
    %end
    %Aex(i,:) = g1.obj';
    g1.modelsense = 'min';
    r = gurobi(g1,params);
    yb(i,1) = r.objval;
    g1.modelsense = 'max';
    r = gurobi(g1,params);
    yb(i,2) = r.objval;
end


% calculating actual confidence intervals
if model.options.ss
    [r,W,J] = stsim(x,model);
else
    [r,W,J] = instsim(x,model);
end
ybactual = yb;
impres = res;
xbest = x;
fbest = impres.fmin;
Aex(:,vrev) = 0;
for i = 1:nconfs
    disp(find(Aex(i,:)));
    Aeq = Aex(i,:)*N1;
    [range,ximp,fimp] = flxlimcalc(x,A,b,actcon,Aeq,model,yb(i,:),r,W,J);
    disp(range);
    ybactual(i,:) = range;
    disp(ybactual);
    if fimp < fbest
        fbest = fimp;
        xbest = ximp;
    end
end
if fbest < impres.fmin
    disp('Improved solution found')
    impres = compileresult(xbest,model);
else
    disp('No improved solution found')
end
% Final FVA for range trimming
%{
c1 = '';
c2 = '';
c1(1:nconfs) = '>';
c2(1:nconfs) = '<';
Aex(:,vrev) = Aex(:,vrev) - Aex(:,vfwd);
gmod1.A = sparse(model.vardata.S_bal);
gmod1.rhs = zeros(size(gmod1.A(:,1)));
gmod1.sense(1:length(gmod1.rhs)) = '=';
gmod1.lb = model.vardata.vb(:,1);
gmod1.ub = model.vardata.vb(:,2);
gmod1.obj = zeros(size(gmod1.lb));
%Aex = zeros(nconfs,length(gmod.lb));
gmod1.vtype(1:length(gmod1.lb)) = 'C';
params.outputflag = 0;
gmod1.A = [gmod1.A;sparse([Aex;Aex])];
gmod1.rhs = [gmod1.rhs;ybactual(:,1);ybactual(:,2)];
gmod1.sense = [gmod1.sense,c1,c2];

for i = 1:length(res.fluxes)
    g1 = gmod1;
    g1.obj(i) = 1;
    if vfwd(i)
        g1.obj(i+1) = -1;
    end
    g1.modelsense = 'min';
    r = gurobi(g1,params);
    res.fluxes(i).vLB = r.objval;
    g1.modelsense = 'max';
    r = gurobi(g1,params);
    res.fluxes(i).vUB = r.objval;
end 

%}

end

