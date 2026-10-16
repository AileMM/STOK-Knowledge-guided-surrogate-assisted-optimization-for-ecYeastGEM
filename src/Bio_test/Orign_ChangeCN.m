load('C_sources.mat');
C_num=length(C_sources);
load('N_sources.mat');
N_num=length(N_sources);
data_pool=[];
i=1;

for j=1:C_num
    for k=1:N_num
        data_pool(i,1)=data(k,j);
        data_pool(i,2)=j;
        data_pool(i,3)=k;
        i=i+1;
    end
end
posGluc = strcmp(ecModel_batch.rxnNames,'D-glucose exchange (reversible)');
posNH4  = strcmp(ecModel_batch.rxnNames,'ammonium exchange (reversible)');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posGluc),0,'u');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posNH4),0,'u');

train_pool=[];
test_pool=[];
chosen=randperm(C_num*N_num);
for i=1:150
    train_pool=[train_pool;data_pool(chosen(i),:)];
end
for i=151:length(chosen)
    test_pool=[test_pool;data_pool(chosen(i),:)];
end

parfor k=1:150
   C_s=train_pool(k,2);
   N_s=train_pool(k,3);
   if strcmp(C_sources{C_s},'(R)-lactate')
       vC = [+Inf,+Inf];
   else
       vC = +Inf;
   end
   ecModel_lim=changeMediaCN(ecModel_batch,C_sources{C_s},N_sources{N_s},vC,+Inf);
   solution=solveLP(ecModel_lim);
end

parfor k=1:150
    C_s=train_pool(k,2);
    N_s=train_pool(k,3);
    if strcmp(C_sources{C_s},'(R)-lactate')
        vC = [+Inf,+Inf];
    else
        vC = +Inf;
    end
    ecModel_lim=changeMediaCN(ecModel_batch,C_sources{C_s},N_sources{N_s},vC,+Inf);
    solution=solveLP(ecModel_lim);
end
    
function model = changeMediaCN(model,C_source,N_source,C_value,N_value)
    model = openUptake(model,C_source,C_value);
    model = openUptake(model,N_source,N_value);
    if strcmp(C_source,'D-fructose')
        model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1134')) = 0;
        model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1134')) = 0;
    elseif strcmp(C_source,'D-mannose')
        model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1139')) = 0;
        model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1139')) = 0;
    end
end

function model = openUptake(model,source,value)
    if strcmp(source,'(R)-lactate')
        model = openUptake(model,'(S)-lactate',value(2));
        value = value(1);
    end
    posB = strcmp(model.rxnNames,[source ' exchange (reversible)']);
    posF = strcmp(model.rxnNames,[source ' exchange']);
    if sum(posB) == 0
        model = changeRxnBounds(model,model.rxns(posF),-value,'l');
    else
        model = changeRxnBounds(model,model.rxns(posB),value,'u');
    end
    model = changeRxnBounds(model,model.rxns(posF),0,'u');
end
               
               
               