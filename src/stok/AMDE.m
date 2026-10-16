clear
clc
load('eny1759.mat');
load('ecYeastGEM_batch.mat');
p=parpool();
p.IdleTimeout=inf;
load('data.mat');
load('C_sources.mat');
C_num=length(C_sources);
load('N_sources.mat');
N_num=length(N_sources);
data_pool=[];
i=1;
% ==========================================================
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
% ==========================================================
protposition=[];
for i=1:length(ecModel_batch.mets)
    if ismember('prot',cell2mat(ecModel_batch.mets(i)))
        protposition=[protposition;i];
    end
end
protposition(length(protposition),:)=[];
for i=1:length(ecModel_batch.rxnNames)
    if ismember('draw_prot_',cell2mat(ecModel_batch.rxnNames(i)))
        drawrxn=i;
        break
    end
end
kcats_inf={};
kcat=[];
k=1;
for i=1:length(protposition)
    posi=find(ecModel_batch.S(protposition(i),:)~=0);
    for j=1:length(posi)
        if posi(j)>drawrxn-1
            continue
        end
        kcats_inf{k}={protposition(i),posi(j),ecModel_batch.S(protposition(i),posi(j))};
        k=k+1;
        kcat=[kcat,ecModel_batch.S(protposition(i),posi(j))];
    end
end
% for j=1:length(kcat)
%     ecModel_batch.S(cell2mat(kcats_inf{j}(1)),cell2mat(kcats_inf{j}(2)))=0;
% end
kcats_inf(eny)=[];
kcat(eny)=[];
kcat=full(kcat);
kcat_lb=kcat*10;
kcat_ub=kcat*0.001;
the_final=[];
best_swarm=[];
error_final=[];
the_best=[];
Times=1;
for time=1:Times
    % ==========================================================
    train_pool=[];
    test_pool=[];
    chosen=randperm(C_num*N_num);
    for i=1:150
        train_pool=[train_pool;data_pool(chosen(i),:)];
    end
    for i=151:length(chosen)
        test_pool=[test_pool;data_pool(chosen(i),:)];
    end
    % ==========================================================
    train_data=train_pool((train_pool(:,1)>0),1);
    test_data=test_pool((test_pool(:,1)>0),1);
    Swarm=[];
    Swarmsize=40;
    for i=1:Swarmsize
        Swarm=[Swarm;kcat_lb+(kcat_ub-kcat_lb).*rand(1,length(kcat))];
    end
    error=[];
    result_swarm=[];
%     Swarm=[Swarm;kcat];
    stra1=0.2*Swarmsize;
    stra2=0.2*Swarmsize;
    stra3=0.2*Swarmsize;
    stra4=0.2*Swarmsize;
    stra5=0.2*Swarmsize;
    iteration=2000;
    for iter=1664:iteration
        tri_swarm=[];
        for i=1:Swarmsize
            if iter==1
                tri_swarm(i,:)=Swarm(randperm(Swarmsize,1),:)+rand(1,length(kcat)).*(Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:));
            elseif (iter<0.2*iteration)&&(iter~=1)
                tri_swarm(i,:)=Swarm(randperm(Swarmsize,1),:)+rand(1,length(kcat)).*(Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:));
            elseif iter>=0.2*iteration
                ran1=randperm(Swarmsize);
                stra_1=ran1(1,1:stra1);
                stra_2=ran1(1,stra1+1:stra1+stra2);
                stra_3=ran1(1,stra1+stra2+1:stra1+stra2+stra3);
                stra_4=ran1(1,stra1+stra2+stra3+1:stra1+stra2+stra3+stra4);
                stra_5=ran1(1,stra1+stra2+stra3+stra4+1:stra1+stra2+stra3+stra4+stra5);
                for u=1:length(stra_1)
                    tri_swarm(stra_1(u),:)=Swarm(randperm(Swarmsize,1),:)+rand(1,length(kcat)).*(Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:));
                end
                for u=1:length(stra_2)
                    tri_swarm(stra_2(u),:)=Swarm(randperm(Swarmsize,1),:)+rand(1,length(kcat)).*(Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:))+rand(1,length(kcat)).*(Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:));
                end
                for u=1:length(stra_3)
                    tri_swarm(stra_3(u),:)=Swarm(row,:)+rand(1,length(kcat)).*(Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:))+rand(1,length(kcat)).*(Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:));
                end
                for u=1:length(stra_4)
                    tri_swarm(stra_4(u),:)=Swarm(stra_4(u),:)+rand(1,length(kcat)).*(Swarm(row,:)-Swarm(stra_4(u),:))+rand(1,length(kcat)).*(Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:));
                end
                for u=1:length(stra_5)
                    tri_swarm(stra_5(u),:)=Swarm(row,:)++rand(1,length(kcat)).*(Swarm(row,:)-Swarm(stra_5(u),:))+rand(1,length(kcat)).*(Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:)+Swarm(randperm(Swarmsize,1),:)-Swarm(randperm(Swarmsize,1),:));
                end
            end
        end
        for i=1:Swarmsize
            for j=1:length(kcat)
                if tri_swarm(i,j)>kcat_ub(j)
                    tri_swarm(i,j)=(kcat_ub(j)-kcat_lb(j))*rand(1,1)+kcat_lb(j);
                elseif tri_swarm(i,j)<kcat_lb(j)
                    tri.swarm(i,j)=(kcat_ub(j)-kcat_lb(j))*rand(1,1)+kcat_lb(j);
                end
            end
        end
        jrand=randi(length(kcat));
        for i=1:Swarmsize
            for j=1:length(kcat)
                if rand(1,1)<0.7||j==jrand
                    child_swarm(i,j)=tri_swarm(i,j);
                else
                    child_swarm(i,j)=Swarm(i,j);
                end
            end
        end
        result_child=[];        
        error_child=[];
        error_swarm=[];
        for i=1:Swarmsize
            for j=1:length(kcat)
                ecModel_batch.S(cell2mat(kcats_inf{j}(1)),cell2mat(kcats_inf{j}(2)))=child_swarm(i,j);
            end
            % ==========================================================
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
                % ==========================================================
               if isempty(solution.f)
                   result_child(i,k)=1;
               else
                   result_child(i,k)=-solution.f;               
               end
            end
            
            error_child(i)=nanmean(abs((train_data(:,1)'-result_child(i,(train_pool(:,1)>0)))./train_data(:,1)'))*100;
            
            if iter==1
                for j=1:length(kcat)
                    ecModel_batch.S(cell2mat(kcats_inf{j}(1)),cell2mat(kcats_inf{j}(2)))=Swarm(i,j);
                end
                % ==========================================================
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
                % ==========================================================
                    if isempty(solution.f)
                        result_swarm(i,k)=1; 
                    else
                        result_swarm(i,k)=-solution.f;
                    end
                end
                error_swarm(i)=nanmean(abs((train_data(:,1)'-result_swarm(i,(train_pool(:,1)>0)))./train_data(:,1)'))*100;
                
                if error_swarm(i)<error_child(i)
                    Swarm(i,:)=Swarm(i,:);
                    error(i)=error_swarm(i);
                else
                    Swarm(i,:)=child_swarm(i,:);
                    error(i)=error_child(i);
                end
            else
                if error(i)<error_child(i)
                    Swarm(i,:)=Swarm(i,:);
                    error(i)=error(i);
                else
                    Swarm(i,:)=child_swarm(i,:);
                    error(i)=error_child(i);
                end
            end
        end   
        
        [m,row]=find(error==min(error));
        error_final(time,iter)=min(error);
        best_swarm(iter,:)=Swarm(row(1),:);
        if iter>=0.2*iteration
            stra_error(1)=mean(error(stra_1));
            stra_error(2)=mean(error(stra_2));
            stra_error(3)=mean(error(stra_3));
            stra_error(4)=mean(error(stra_4));
            stra_error(5)=mean(error(stra_5));
            [~,row_min]=find(stra_error==min(stra_error));
            [~,row_max]=find(stra_error==max(stra_error));
            if row_max==1
                if stra1==1
                    if stra2==max([stra1,stra2,stra3,stra4,stra5])
                        stra2=stra2-1;
                    elseif stra3==max([stra1,stra2,stra3,stra4,stra5])
                        stra3=stra3-1;
                    elseif stra4==max([stra1,stra2,stra3,stra4,stra5])
                        stra4=stra4-1;
                    elseif stra5==max([stra1,stra2,stra3,stra4,stra5])
                        stra5=stra5-1;
                    end
                    if row_min==2
                        stra2=stra2+1;
                    elseif row_min==3
                        stra3=stra3+1;
                    elseif row_min==4
                        stra4=stra4+1;
                    elseif row_min==5
                        stra5=stra5+1;
                    end
                else
                    stra1=stra1-1;
                    if row_min==2
                        stra2=stra2+1;
                    elseif row_min==3
                        stra3=stra3+1;
                    elseif row_min==4
                        stra4=stra4+1;
                    elseif row_min==5
                        stra5=stra5+1;
                    end
                end
            elseif row_max==2
                if stra2==1
                    if stra1==max([stra1,stra2,stra3,stra4,stra5])
                        stra1=stra1-1;
                    elseif stra3==max([stra1,stra2,stra3,stra4,stra5])
                        stra3=stra3-1;
                    elseif stra4==max([stra1,stra2,stra3,stra4,stra5])
                        stra4=stra4-1;
                    elseif stra5==max([stra1,stra2,stra3,stra4,stra5])
                        stra5=stra5-1;
                    end
                    if row_min==1
                        stra1=stra1+1;
                    elseif row_min==3
                        stra3=stra3+1;
                    elseif row_min==4
                        stra4=stra4+1;
                    elseif row_min==5
                        stra5=stra5+1;
                    end
                else
                    stra2=stra2-1;
                    if row_min==1
                        stra1=stra1+1;
                    elseif row_min==3
                        stra3=stra3+1;
                    elseif row_min==4
                        stra4=stra4+1;
                    elseif row_min==5
                        stra5=stra5+1;
                    end
                end
            elseif row_max==3
                if stra3==1
                    if stra2==max([stra1,stra2,stra3,stra4,stra5])
                        stra2=stra2-1;
                    elseif stra1==max([stra1,stra2,stra3,stra4,stra5])
                        stra1=stra1-1;
                    elseif stra4==max([stra1,stra2,stra3,stra4,stra5])
                        stra4=stra4-1;
                    elseif stra5==max([stra1,stra2,stra3,stra4,stra5])
                        stra5=stra5-1;
                    end
                    if row_min==2
                        stra2=stra2+1;
                    elseif row_min==1
                        stra1=stra1+1;
                    elseif row_min==4
                        stra4=stra4+1;
                    elseif row_min==5
                        stra5=stra5+1;
                    end
                else
                    stra3=stra3-1;
                    if row_min==2
                        stra2=stra2+1;
                    elseif row_min==1
                        stra1=stra1+1;
                    elseif row_min==4
                        stra4=stra4+1;
                    elseif row_min==5
                        stra5=stra5+1;
                    end
                end
            elseif row_max==4
                if stra4==1
                    if stra2==max([stra1,stra2,stra3,stra4,stra5])
                        stra2=stra2-1;
                    elseif stra3==max([stra1,stra2,stra3,stra4,stra5])
                        stra3=stra3-1;
                    elseif stra1==max([stra1,stra2,stra3,stra4,stra5])
                        stra1=stra1-1;
                    elseif stra5==max([stra1,stra2,stra3,stra4,stra5])
                        stra5=stra5-1;
                    end
                    if row_min==2
                        stra2=stra2+1;
                    elseif row_min==3
                        stra3=stra3+1;
                    elseif row_min==1
                        stra1=stra1+1;
                    elseif row_min==5
                        stra5=stra5+1;
                    end
                else
                    stra4=stra4-1;
                    if row_min==2
                        stra2=stra2+1;
                    elseif row_min==3
                        stra3=stra3+1;
                    elseif row_min==1
                        stra1=stra1+1;
                    elseif row_min==5
                        stra5=stra5+1;
                    end
                end
            elseif row_max==5
                if stra5==1
                    if stra2==max([stra1,stra2,stra3,stra4,stra5])
                        stra2=stra2-1;
                    elseif stra3==max([stra1,stra2,stra3,stra4,stra5])
                        stra3=stra3-1;
                    elseif stra4==max([stra1,stra2,stra3,stra4,stra5])
                        stra4=stra4-1;
                    elseif stra1==max([stra1,stra2,stra3,stra4,stra5])
                        stra1=stra1-1;
                    end
                    if row_min==2
                        stra2=stra2+1;
                    elseif row_min==3
                        stra3=stra3+1;
                    elseif row_min==4
                        stra4=stra4+1;
                    elseif row_min==1
                        stra1=stra1+1;
                    end
                else
                    stra5=stra5-1;
                    if row_min==2
                        stra2=stra2+1;
                    elseif row_min==3
                        stra3=stra3+1;
                    elseif row_min==4
                        stra4=stra4+1;
                    elseif row_min==1
                        stra1=stra1+1;
                    end
                end
            end
        end
        save 'result2502_4.mat';
    end
    plot(1:iteration,error_final(time,:));
    hold on
    [final,row]=find(error_final(time,:)==min(error_final(time,:)));
    the_best(time,:)=best_swarm(row(1),:);
    the_final(time)=min(error_final(time,:));
 end

delete(p);

% for j=1:length(kcat)
%     ecModel_batch.S(cell2mat(kcats_inf{j}(1)),cell2mat(kcats_inf{j}(2)))=the_best(time,j);
% end
% save('result.mat','ecModel_batch');


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
