%% 
citys= [[1304,2312];[3639,1315];[4177,2244];[3712,1399];[3488,1535];[3326,1556];[3238,1229];[4196,1004];[4312,790];[4386,570];[3007,1970];[2562,1756];[2788,1491];[2381,1676];[1332,695];[3715,1678];[3918,2179];[4061,2370];[3780,2212];[3676,2578];[4029,2838];[4263,2931];[3429,1908];[3507,2367];[3394,2643];[3439,3201];[2935,3240];[3140,3550];[2545,2357];[2778,2826];[2370,2975]];

n = size(citys,1);% 31个城市
%% 计算城市之间的距离
D = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            D(i,j) = sqrt(((citys(i,1) - citys(j,1))^2)+((citys(i,2) - citys(j,2))^2));
        else
            D(i,j) = eps;
        end
    end
end

m = 50;                              % 蚂蚁数量
alpha = 1;                           % 信息素重要程度因子
beta = 5;                            % 启发函数重要程度因子
rho = 0.1;                           % 信息素挥发因子
Q = 1;                               % 常系数（信息素释放量）
Eta = 1./D;                          % 启发函数   // 距离作为启发函数
Tau = ones(n,n);                     % 信息素矩阵
Table = zeros(m,n);                  % 路径记录表 50x31的矩阵
iter = 1;                            % 迭代次数初值
iter_max = 200;                      % 最大迭代次数
Route_best = zeros(iter_max,n);      % 各代最佳路径
Length_best = zeros(iter_max,1);     % 各代最佳路径的长度
Length_ave = zeros(iter_max,1);      % 各代路径的平均长度


while iter <= iter_max
    % 随机产生各个蚂蚁的起点城市  m是蚂蚁的数量
      start = zeros(m,1);
      for i = 1:m   
          temp = randperm(n);  % n是城市的个数  随机选取
          start(i) = temp(1);
      end
      Table(:,1) = start; %出发地作为第1列  51x30矩阵
      % 构建解空间
      citys_index = 1:n;
      % 逐个蚂蚁路径选择
      for i = 1:m
          % 逐个城市路径选择
         for j = 2:n
             tabu = Table(i,1:(j - 1)); % 已访问的城市集合 (1:j) 相当于subString
             allow_index = ~ismember(citys_index,tabu); %取反 0-1   1-0
             allow = citys_index(allow_index);% 待访问的城市集合
             P = allow;
             % 计算城市间转移概率
             for k = 1:length(allow)
                 %%  Tau 信息素  Eta启发函数（以2城市之间的欧式距离的倒数） end是最后一列
                 %%  tabu(end)刚访问过的城市编号  tau 记录已访问的城市集合（46行，代码风格不好，起名visitedCity）
                 P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
             end
             P = P/sum(P);
             % 轮盘赌法选择下一个访问城市
             Pc = cumsum(P);%[2:n] 依次累加
            target_index = find(Pc >= rand); %大于阈值随机选
            target = allow(target_index(1));
            Table(i,j) = target;
         end
      end
      % 计算各个蚂蚁的路径距离
      Length = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);%% 蚂蚁走的路径对应的城市编号  i,表示取第i行  ：表示取所有列
          for j = 1:(n - 1)
              Length(i) = Length(i) + D(Route(j),Route(j + 1));
          end
          %% 最后一个节点到原点单独计算
          Length(i) = Length(i) + D(Route(n),Route(1));
      end
      % 计算最短路径距离及平均距离
      if iter == 1
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min_Length;
          Length_ave(iter) = mean(Length);
          Route_best(iter,:) = Table(min_index,:);
      else
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min(Length_best(iter - 1),min_Length);%上一轮迭代最优与本轮值取较小
          Length_ave(iter) = mean(Length);
          if Length_best(iter) == min_Length
              Route_best(iter,:) = Table(min_index,:);%本轮迭代为最优路径
          else
              Route_best(iter,:) = Route_best((iter-1),:);% 保留上一轮最优路径作为本轮的最优路径
          end
      end
      % 更新信息素
      Delta_Tau = zeros(n,n);
      % 逐个蚂蚁计算
      for i = 1:m
          % 逐个城市计算
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
      end
      Tau = (1-rho) * Tau + Delta_Tau;
    % 迭代次数加1，清空路径记录表
    iter = iter + 1;
    Table = zeros(m,n);
end


[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['最短距离:' num2str(Shortest_Length)]);
disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);


figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       起点');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       终点');
xlabel('城市位置横坐标')
ylabel('城市位置纵坐标')
title(['蚁群算法优化路径(最短距离:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r:')
legend('最短距离','平均距离')
xlabel('迭代次数')
ylabel('距离')
title('各代最短距离与平均距离对比')
