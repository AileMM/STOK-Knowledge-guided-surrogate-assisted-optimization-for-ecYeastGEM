function problem= problem_family(problem_id,M,d)

switch(problem_id)
    case 1 
        dim = d;
        problem.M = M;
        problem.name = 'DTLZ1';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = UniformPoint(500,M)/2;
        problem.fnc = @(x)DTLZ1(x,problem.M,optimum);
        problem.F_best = optimum;
    case 2 
        dim = d;
        problem.M = M;
        problem.name = 'DTLZ2';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        
        optimum = UniformPoint(500,M);
        optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,M);
        
        problem.fnc = @(x)DTLZ2(x,problem.M,optimum);
        problem.F_best = optimum;
    case 3 
        dim = d;
        problem.M = M;
        problem.name = 'DTLZ3';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = UniformPoint(500,M);
        optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,M);
        problem.fnc = @(x)DTLZ3(x,problem.M,optimum);
        problem.F_best = optimum;
    case 4 
        dim = d;
        problem.M = M;
        problem.name = 'DTLZ4';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = UniformPoint(500,M);
        optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,M);
        problem.fnc = @(x)DTLZ4(x,problem.M,optimum);
        problem.F_best = optimum;
    case 5 
        dim = d;
        problem.M = M;
        problem.name = 'DTLZ5';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        
        optimum = [0:1/(500-1):1;1:-1/(500-1):0]';
        optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,size(optimum,2));
        optimum = [optimum(:,ones(1,M-2)),optimum];
        optimum = optimum./sqrt(2).^repmat([M-2,M-2:-1:0],size(optimum,1),1);
        rotation = eye(dim,dim);
        problem.fnc = @(x)DTLZ5(x,problem.M,optimum);
        problem.F_best = optimum;
    case 6 
        dim = d;
        problem.M = M;
        problem.name = 'DTLZ6';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = [0:1/(500-1):1;1:-1/(500-1):0]';
        optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,size(optimum,2));
        optimum = [optimum(:,ones(1,M-2)),optimum];
        optimum = optimum./sqrt(2).^repmat([M-2,M-2:-1:0],size(optimum,1),1);
        rotation = eye(dim,dim);
        problem.fnc = @(x)DTLZ6(x,problem.M,optimum);
        problem.F_best = optimum;
    case 7 
        dim = d;
        problem.M = M;
        problem.name = 'DTLZ7';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        
        interval     = [0,0.251412,0.631627,0.859401];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
        X            = UniformPoint(500,M-1,'grid');
        X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        optimum            = [X,2*(M-sum(X/2.*(1+sin(3*pi.*X)),2))];
        rotation = eye(dim,dim);
        problem.fnc = @(x)DTLZ7(x,problem.M,optimum);
        problem.F_best = optimum;
    case 8 
        dim = d;
        problem.M = M;
        problem.name = 'WFG1';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = W1_Optimum(500,M);
        problem.fnc = @(x)WFG1(x,problem.M,optimum);
        problem.F_best = optimum;
    case 9 
        dim = d;
        problem.M = M;
        problem.name = 'WFG2';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = W2_Optimum(500,M);
        problem.fnc = @(x)WFG2(x,problem.M,optimum);
        problem.F_best = optimum;
    case 10 
        dim = d;
        problem.M = M;
        problem.name = 'WFG3';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = W3_Optimum(500,M);
        problem.fnc = @(x)WFG3(x,problem.M,optimum);
        problem.F_best = optimum;
    case 11
        dim = d;
        problem.M = M;
        problem.name = 'WFG4';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = W4_Optimum(500,M);
        problem.fnc = @(x)WFG4(x,problem.M,optimum);
        problem.F_best = optimum;
    case 12 
        dim = d;
        problem.M = M;
        problem.name = 'WFG5';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = W5_Optimum(500,M);
        problem.fnc = @(x)WFG5(x,problem.M,optimum);
        problem.F_best = optimum;
    case 13 
        dim = d;
        problem.M = M;
        problem.name = 'ZDT1';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = Z1_Optimum(500,M);
        problem.fnc = @(x)ZDT1(x,problem.M,optimum);
        problem.F_best = optimum;
    case 14 
        dim = d;
        problem.M = M;
        problem.name = 'ZDT2';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = Z2_Optimum(500,M);
        problem.fnc = @(x)ZDT2(x,problem.M,optimum);
        problem.F_best = optimum;
    case 15 
        dim = d;
        problem.M = M;
        problem.name = 'ZDT3';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = Z3_Optimum(500,M);
        problem.fnc = @(x)ZDT3(x,problem.M,optimum);
        problem.F_best = optimum;
    case 16 
        dim = d;
        problem.M = M;
        problem.name = 'ZDT4';
        problem.dims = dim;
        problem.lb = [0, -5*ones(1, dim-1)];
        problem.ub = [1,  5*ones(1, dim-1)];
        optimum = Z4_Optimum(500,M);
        problem.fnc = @(x)ZDT4(x,problem.M,optimum);
        problem.F_best = optimum;
    case 17 
        dim = d;
        problem.M = M;
        problem.name = 'ZDT6';
        problem.dims = dim;
        problem.lb = 0*ones(1,dim);
        problem.ub = 1*ones(1,dim);
        optimum = Z6_Optimum(500,M);
        problem.fnc = @(x)ZDT6(x,problem.M,optimum);
        problem.F_best = optimum;
   
end
