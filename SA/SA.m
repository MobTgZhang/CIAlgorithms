clear;
clc;
%% preparing data set
filename = 'kroA200.tsp';
url = 'http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/kroA200.tsp';
if ~exist(filename)
    urlwrite(url,filename);
end
file = fopen(filename);
flag = false;
citys_mat = [];
while ~feof(file)
  str = fgetl(file);
  if strfind(str,'EOF')
    break
  elseif strfind(str,'NODE_COORD_SECTION')
    flag = true;
  elseif flag
    a = strsplit(str);
    index = str2num(char(a(1)));
    in_x = str2num(char(a(2)));
    in_y = str2num(char(a(3)));
    if ~citys_mat
      citys_mat = [in_x,in_y];
    end
    citys_mat = [citys_mat;[in_x,in_y]];
  else
    continue
  end
end
%% 将城市打印在图上
x = citys_mat(:,1);
y = citys_mat(:,2);
% figure;
% set(gcf,'position',[20,50,1500,800]);
% scatter(x,y,'r');
% grid();
%% 给定参数值
alpha = 0.80;
t0 = 200;
tf = 1;
Markov_length = 10000;
%% 计算城市之间的距离矩阵
citys_length = length(citys_mat);
citys_x = citys_mat(:,1) * ones(1,citys_length);
citys_y = citys_mat(:,2) * ones(1,citys_length);
citys_distance_mat = sqrt((citys_x-citys_x').^2+(citys_y-citys_y').^2);
sol_new = 1:citys_length;
sol_best = sol_new;
E_current = inf;
E_best = inf;
sol_current = sol_new;sol_beat = sol_new;
t = t0;
epoches = 0;
length_list = [];
while t>tf
    %产生随机扰动
    for k=1:Markov_length
        if(rand<0.5)
            %两个值进行交换
            ind1=0;ind2=0;
            while (ind1==ind2)
                ind1 = ceil(rand .* citys_length);
                ind2 = ceil(rand .* citys_length);
            end
            tmp = sol_new(ind1);
            sol_new(ind1) = sol_new(ind2);
            sol_new(ind2) = tmp;
        else
            %三交换
            ind1=0;ind2=0;ind3=0;
            while(ind1==ind2)||(ind1==ind3)||(ind2==ind3)
                ind1 = ceil(rand .* citys_length);
                ind2 = ceil(rand .* citys_length);
                ind3 = ceil(rand .* citys_length);
            end
            tmp1=ind1;tmp2 = ind2;tmp3=ind3;
            %确保ind1<ind2<ind3
            if((ind1<ind2)&&(ind2<ind3))
            elseif (ind1<ind3)&&(ind3<ind2)
                ind2 = tmp3;ind3=tmp2;
            elseif (ind2<ind1)&&(ind1<ind3)
                ind1 = tmp2;ind2 = tmp1;ind3 = tmp3;
            elseif (ind2<ind3)&&(ind3<ind1)
                ind1 = tmp2;ind2 = tmp3;ind3 = tmp1;
            elseif (ind3<ind1)&&(ind1<ind2)
                ind1 = tmp3;ind2 = tmp1;ind3 = tmp2;
            elseif (ind3<ind2)&&(ind2<ind1)
                ind1 = tmp3;ind2 = tmp2;ind3 = tmp1;
            end
            tmplist = sol_new((ind1+1):(ind2-1));
            sol_new((ind1+1):(ind1+ind3-ind2+1)) = sol_new(ind2:ind3);
            sol_new((ind1+ind3-ind2+2):ind3) = tmplist;
        end
        E_new = 0;
        for i=1:(citys_length-1)
            E_new = E_new +citys_distance_mat(sol_new(i),sol_new(i+1));
        end
        E_new = E_new +citys_distance_mat(sol_new(citys_length),sol_new(1));
        if E_new < E_current
            E_current = E_new;
            sol_current = sol_new;
            if E_new < E_best
                E_best = E_new;
                sol_best = sol_new;
            end
        else
            %若新的解目标函数小于当前的解，则以一定的概率接受新解
            if rand <exp(-(E_new-E_current)./t)
                E_current = E_new;
                sol_current = sol_new;
            else
                sol_new = sol_current;
            end
        end
    end
    length_list(end+1) = E_best;
    t = t*alpha;
    epoches = epoches +1;
    fprintf('Epoches:%d,temperature:%.2f\n',epoches,t);
end
for k=1:citys_length-1
    dotA = citys_mat(sol_best(k),:);
    dotB = citys_mat(sol_best(k+1),:);
    scatter(dotA(1),dotA(2),'b');
    scatter(dotB(1),dotB(2),'b');
    plot([dotA(1),dotB(1)],[dotA(2),dotB(2)],'color','r');
    hold on;
end
dotA = citys_mat(citys_length,:);
dotB = citys_mat(1,:);
scatter(dotA(1),dotA(2),'b');
scatter(dotB(1),dotB(2),'b');
plot([dotA(1),dotB(1)],[dotA(2),dotB(2)],'color','r');
saveas(gcf,'route.jpg');
close(gcf);
x = 1:length(length_list);
plot(x,length_list);
saveas(gcf,'loss.jpg');

