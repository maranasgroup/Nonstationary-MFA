function [res] = flxestimate(model,res)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = model.vardata.N;
nu = model.vardata.nu;
if ~model.options.ss
    nc = model.vardata.npool;
else
    nc = 0;
end
nh = sum([model.data.nh]);
A = [N;-N];
A = blkdiag(A,eye(nc+nh));
b = model.vardata.vb;
b(:,2) = -1*b(:,2);
b = b(:);
b = [b;1e-5*ones(nc,1);1e-7*ones(nh,1)];
nms = model.options.multistart;
fopt = Inf;
for i = 1:nms
    if model.options.reinit
        [x,actcon] = initialize(model);
    else
        if model.options.ss
            x = [res.reinit_data.u;res.reinit_data.h];
        else
            x = [res.reinit_data.u;res.reinit_data.c;res.reinit_data.h];
        end
        if any(A*x < b)
            options = optimset('Display','Iter','MaxFunEvals',100000,'maxiter',10000,'Algorithm','interior-point');
            u0 = x(1:nu);
            v = N*u0;
            A1 = [-N;N];
            b1 = [-model.vardata.vb(:,1);model.vardata.vb(:,2)];
            lb1 = 1e-7*ones(nu,1);
            ub1 = 1e8*ones(nu,1);
            fun1=@(u)(v-(N*u))'*(v-(N*u));
            [u,~,exitflag]=fmincon(fun1,u0,A1,b1,[],[],lb1,ub1,[],options);
            x(1:nu) = u;
        end
        actcon = find(A*x<=b);
    end
    [xf,f,~,actcon] = lsqsolve(x,model,A,b,actcon);
    if f<fopt
        fopt = f;
        xopt = xf;
    end
end

model.options.dfbase = eps;
iter = 1;
fail = true;
x = xopt;
f = fopt;
ac = actcon;
while fail || iter <= 5
    [x,f,fail,ac] = lsqsolve(x,model,A,b,ac);
    if f<fopt
        xopt = x;
        fopt = f;
        actcon = ac;
    end
    iter = iter+1;
end
res = compileresult(xopt,model);
end

