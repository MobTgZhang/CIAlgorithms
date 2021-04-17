
import numpy as np
import matplotlib.pyplot as plt
import urllib.request
import os
import time
def download(root_path,filename):
    if not os.path.exists(root_path):
        os.mkdir(root_path)
    if not os.path.exists(os.path.join(root_path,filename)):
        url = "http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/"+filename
        urllib.request.urlretrieve(url,os.path.join(root_path,filename))
        print("The data set: %s downloaded!"%os.path.join(root_path,filename))
    else:
        print("The data set: %s already has downloaded!"%os.path.join(root_path,filename))
def get_data(filename):
    data_list = []
    with open(filename,mode="r") as f:
        flag = False
        while True:
            line = f.readline()
            if "EOF" in line:
                break
            elif "NODE_COORD_SECTION" in line:
                flag = True
            elif flag:
                tmp = line.strip().split(" ")
                data_list.append([float(item) for item in tmp])
    return np.array(data_list)
class ACO:
    def __init__(self,ant_num,alpha,beta,rho,Q,epoches):
        self.ant_num = ant_num
        self.alpha = alpha
        self.beta = beta
        self.rho = rho
        self.Q = Q
        self.epoches = epoches
        self.citys_mat = None
        self.E_best = None
        self.sol_best = None
        self.length_list = None
        self.name = time.strftime("%Y%m%d%H%M", time.localtime(time.time()))
    def solve(self,citys_mat):
        self.citys_mat = citys_mat
        citys_num = citys_mat.shape[0]
        # 获取邻接矩阵
        citys_x = citys_mat[:, 0].reshape(citys_num, 1).dot(np.ones((1, citys_num)))
        citys_y = citys_mat[:, 1].reshape(citys_num, 1).dot(np.ones((1, citys_num)))
        citys_distance = np.sqrt(np.square(citys_x - citys_x.T) + np.square(citys_y - citys_y.T))
        # 初始化启发函数
        Heu_f = 1.0/(citys_distance + np.diag([np.inf] * citys_num))
        # 信息素矩阵
        Tau_table = np.ones((citys_num,citys_num))
        # 每一次迭代过程中每个蚂蚁的路径记录表
        Route_table = np.zeros((self.ant_num,citys_num),dtype=np.int)
        # 每一次迭代过程中的最佳路径
        Route_best = np.zeros((self.epoches,citys_num),dtype=np.int)
        # 每一次迭代过程中最佳路径记录表
        Length_best = np.zeros(self.epoches)
        # 每次迭代过程中蚂蚁的平均路径长度
        Length_average = np.zeros(self.epoches)
        # 每次迭代过程中当前路径长度
        Length_current = np.zeros(self.ant_num)
        iter = 0
        while iter <self.epoches:
            # 产生城市集合表
            # 随机产生各个蚂蚁的起点城市
            Route_table[:,0]= self.randseed(citys_num)
            # 更新信息素
            Delta_tau = np.zeros((citys_num, citys_num))
            for k in range(self.ant_num):
                # 用于记录蚂蚁下一个访问的城市集合
                # 蚂蚁已经访问过的城市
                tabu = [Route_table[k,0]]
                allow_set = list(set(range(citys_num))-set(tabu))
                city_index = Route_table[k,0]
                for i in range(1,citys_num):
                    # 初始化城市之间的转移概率
                    P_table = np.zeros(len(allow_set))
                    # 计算城市之间的转移概率
                    for j in range(len(allow_set)):
                        P_table[j] = np.power(Tau_table[city_index,allow_set[j]],self.alpha)*\
                                     np.power(Heu_f[city_index,allow_set[j]],self.beta)
                    P_table = P_table/np.sum(P_table)

                    # 轮盘赌算法来选择下一个访问的城市
                    #out_prob = np.cumsum(P_table)
                    while True:
                        r = np.random.rand()
                        index_need = np.where(P_table > r)[0]
                        if len(index_need) >0:
                            city_index2 = allow_set[index_need[0]]
                            break
                    Route_table[k,i]  = city_index2
                    tabu.append(city_index2)
                    allow_set = list(set(range(0,citys_num))-set(tabu))
                    city_index = city_index2
                tabu.append(tabu[0])
                # 计算蚂蚁路径的距离信息
                for j in range(citys_num):
                    Length_current[k] = Length_current[k] + citys_distance[tabu[j],tabu[j+1]]
                for j in range(citys_num):
                    Delta_tau[tabu[j],tabu[j+1]] = Delta_tau[tabu[j],tabu[j+1]] + self.Q / Length_current[k]
            # 计算最短路径、最短路径长度以及平均路径长度
            Length_best[iter] = np.min(Length_current)
            index = np.where(Length_current == np.min(Length_current))[0][0]
            Route_best[iter] = Route_table[index]
            Length_average[iter] = np.mean(Length_current)
            #更新信息素
            Tau_table = (1-self.rho)*Tau_table + Delta_tau
            #Route_table = np.zeros((self.ant_num,citys_num),dtype=np.int)
            Length_current = np.zeros(self.ant_num)

            print("epoches:%d,best value every epoches%.4f"%(iter, Length_best[iter]))
            iter = iter + 1
        self.E_best = np.min(Length_best)
        index = np.where(Length_best == np.min(Length_best))[0][0]
        self.sol_best = Route_table[index]
        self.length_list = Length_average
    def randseed(self,citys_num):
        if self.ant_num <citys_num:
            initial_route = np.random.permutation(range(citys_num))[:self.ant_num]
        else:
            initial_route = np.zeros((self.ant_num,))
            initial_route[:citys_num] = np.random.permutation(range(citys_num))
            tmp_index = citys_num
            while tmp_index + citys_num <= self.ant_num:
                initial_route[tmp_index:citys_num + tmp_index] = np.random.permutation(range(citys_num))
                tmp_index += citys_num
            tmp_left = self.ant_num % citys_num
            if tmp_left != 0:
                initial_route[tmp_index:] = np.random.permutation(range(citys_num))[:tmp_left]
        return initial_route
    def draw(self):
        print(self.sol_best)
        print(self.E_best)
        if not os.path.exists("log"):
            os.mkdir("log")
        # draw loss
        x = np.linspace(0, len(self.length_list) - 1, len(self.length_list))
        y = np.array(self.length_list)
        plt.plot(x, y)
        plt.title(label="loss")
        plt.savefig(os.path.join("log", "%s_loss.png" % self.name))
        plt.close()
        # draw dots
        for k in range(0, len(self.sol_best) - 1):
            start = self.citys_mat[self.sol_best[k]]
            end = self.citys_mat[self.sol_best[k + 1]]
            plt.plot(start[0], start[1], "bo")
            plt.plot(end[0], end[1], "bo")
            plt.arrow(start[0], start[1], end[0] - start[0], end[1] - start[1],
                      length_includes_head=True, head_width=0.2, head_length=0.3, lw=1,
                      color="r")
        start = self.citys_mat[self.sol_best[-1]]
        end = self.citys_mat[self.sol_best[0]]
        plt.plot(start[0], start[1], "bo")
        plt.plot(end[0], end[1], "bo")
        plt.arrow(start[0], start[1], end[0] - start[0], end[1] - start[1],
                  length_includes_head=True, head_width=0.2, head_length=0.3, lw=1,
                  color="r")
        plt.title(label="length:%.2f" % (self.E_best))
        plt.savefig(os.path.join("log", "%s_route.png" % self.name))
        plt.show()
def main():
    filename = "eil51.tsp"
    root_path = "data"
    download(root_path,filename)
    data_list = get_data(os.path.join(root_path,filename))
    ant_num = 500
    alpha = 1
    beta = 5
    rho = 0.2
    Q = 10
    epoches = 20
    model = ACO(ant_num, alpha, beta, rho, Q, epoches)
    model.solve(data_list[:,1:])
    model.draw()
if __name__ == '__main__':
    main()

