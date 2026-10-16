function [solution_externel,imp_max] = knowledge_transfer(database,lb,ub,knowledge_base,surrogates_source,ada_vectors)

num_sources = length(knowledge_base); % the number of solved source tasks
dim = length(lb);
solutions_target = (database(:,1:end-1)-repmat(lb,size(database,1),1))./(repmat(ub,size(database,1),1)-repmat(lb,size(database,1),1));
objs_target = database(:,end);
solutions = zeros(num_sources,dim);
improvements = zeros(num_sources,1);
similarities = zeros(num_sources,1);
ranks_target = zeros(length(objs_target),1);
for i = 1:length(objs_target)  %database中的每个适应度值在其中的排名
    ranks_target(i) = sum(objs_target<objs_target(i))+1;
end

for i = 1:num_sources
    solutions_source = knowledge_base(i).solutions;% 生成的5个源任务的个体是不一致的
    objs_source = knowledge_base(i).objs;
    %这里是每一个源任务的代理模型都去评估目标任务，但是我们的思想是先找到一个最相似的任务直接评估？
    objs_val = zeros(length(objs_target),1); 
    %源任务代理模型预测目标任务解表现，使用源任务的代理模型预测目标任务解的目标值
    for j = 1:length(objs_target)  
        objs_val(j) = surrogates_source(i).func(solutions_target(j,:)-ada_vectors(i,:));
    end
    ranks_val = zeros(length(objs_target),1);
    for j = 1:length(objs_target)
        ranks_val(j) = sum(objs_val<objs_val(j))+1; % 计算基于源任务代理模型预测的objs_val 在目标任务中的排名
    end
    uncertainty = 1-max(abs(ada_vectors(i,:))); % 计算源任务代理模型预测的排名和真实目标任务排名的斯皮尔曼相关系数，衡量源任务代理模型对目标任务的适用性
    similarities(i) = uncertainty*mySpearman(ranks_target,ranks_val);
    %判断任务库中哪些源任务对目标任务有更高的迁移价值，后面还需映射
    improvements(i) = improvement_external(sort(objs_source,'descend'),sort(objs_val,'descend'),sort(objs_target,'descend'),similarities(i));
    [~,idx] = min(objs_source);
    x_ada_source = (solutions_source(idx,:)+ada_vectors(i,:));
    x_ada_source(x_ada_source<0) = 0;
    x_ada_source(x_ada_source>1) = 1;
    solutions(i,:) = lb+x_ada_source.*(ub-lb);
end
[imp_max,idx] = max(improvements);
solution_externel = solutions(idx,:);
% fprintf('/*Solution transfer under the highest rank similarity with %.2f [%.2f to %.2f]*/\n',max(similarities), min(similarities), max(similarities));