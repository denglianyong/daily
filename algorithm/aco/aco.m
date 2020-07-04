%% 
citys= [[1304,2312];[3639,1315];[4177,2244];[3712,1399];[3488,1535];[3326,1556];[3238,1229];[4196,1004];[4312,790];[4386,570];[3007,1970];[2562,1756];[2788,1491];[2381,1676];[1332,695];[3715,1678];[3918,2179];[4061,2370];[3780,2212];[3676,2578];[4029,2838];[4263,2931];[3429,1908];[3507,2367];[3394,2643];[3439,3201];[2935,3240];[3140,3550];[2545,2357];[2778,2826];[2370,2975]];

n = size(citys,1);% 31������
%% �������֮��ľ���
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

m = 50;                              % ��������
alpha = 1;                           % ��Ϣ����Ҫ�̶�����
beta = 5;                            % ����������Ҫ�̶�����
rho = 0.1;                           % ��Ϣ�ػӷ�����
Q = 1;                               % ��ϵ������Ϣ���ͷ�����
Eta = 1./D;                          % ��������   // ������Ϊ��������
Tau = ones(n,n);                     % ��Ϣ�ؾ���
Table = zeros(m,n);                  % ·����¼�� 50x31�ľ���
iter = 1;                            % ����������ֵ
iter_max = 200;                      % ����������
Route_best = zeros(iter_max,n);      % �������·��
Length_best = zeros(iter_max,1);     % �������·���ĳ���
Length_ave = zeros(iter_max,1);      % ����·����ƽ������


while iter <= iter_max
    % ��������������ϵ�������  m�����ϵ�����
      start = zeros(m,1);
      for i = 1:m   
          temp = randperm(n);  % n�ǳ��еĸ���  ���ѡȡ
          start(i) = temp(1);
      end
      Table(:,1) = start; %��������Ϊ��1��  51x30����
      % ������ռ�
      citys_index = 1:n;
      % �������·��ѡ��
      for i = 1:m
          % �������·��ѡ��
         for j = 2:n
             tabu = Table(i,1:(j - 1)); % �ѷ��ʵĳ��м��� (1:j) �൱��subString
             allow_index = ~ismember(citys_index,tabu); %ȡ�� 0-1   1-0
             allow = citys_index(allow_index);% �����ʵĳ��м���
             P = allow;
             % ������м�ת�Ƹ���
             for k = 1:length(allow)
                 %%  Tau ��Ϣ��  Eta������������2����֮���ŷʽ����ĵ����� end�����һ��
                 %%  tabu(end)�շ��ʹ��ĳ��б��  tau ��¼�ѷ��ʵĳ��м��ϣ�46�У������񲻺ã�����visitedCity��
                 P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
             end
             P = P/sum(P);
             % ���̶ķ�ѡ����һ�����ʳ���
             Pc = cumsum(P);%[2:n] �����ۼ�
            target_index = find(Pc >= rand); %������ֵ���ѡ
            target = allow(target_index(1));
            Table(i,j) = target;
         end
      end
      % ����������ϵ�·������
      Length = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);%% �����ߵ�·����Ӧ�ĳ��б��  i,��ʾȡ��i��  ����ʾȡ������
          for j = 1:(n - 1)
              Length(i) = Length(i) + D(Route(j),Route(j + 1));
          end
          %% ���һ���ڵ㵽ԭ�㵥������
          Length(i) = Length(i) + D(Route(n),Route(1));
      end
      % �������·�����뼰ƽ������
      if iter == 1
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min_Length;
          Length_ave(iter) = mean(Length);
          Route_best(iter,:) = Table(min_index,:);
      else
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min(Length_best(iter - 1),min_Length);%��һ�ֵ��������뱾��ֵȡ��С
          Length_ave(iter) = mean(Length);
          if Length_best(iter) == min_Length
              Route_best(iter,:) = Table(min_index,:);%���ֵ���Ϊ����·��
          else
              Route_best(iter,:) = Route_best((iter-1),:);% ������һ������·����Ϊ���ֵ�����·��
          end
      end
      % ������Ϣ��
      Delta_Tau = zeros(n,n);
      % ������ϼ���
      for i = 1:m
          % ������м���
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
      end
      Tau = (1-rho) * Tau + Delta_Tau;
    % ����������1�����·����¼��
    iter = iter + 1;
    Table = zeros(m,n);
end


[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['��̾���:' num2str(Shortest_Length)]);
disp(['���·��:' num2str([Shortest_Route Shortest_Route(1)])]);


figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       ���');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       �յ�');
xlabel('����λ�ú�����')
ylabel('����λ��������')
title(['��Ⱥ�㷨�Ż�·��(��̾���:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r:')
legend('��̾���','ƽ������')
xlabel('��������')
ylabel('����')
title('������̾�����ƽ������Ա�')
