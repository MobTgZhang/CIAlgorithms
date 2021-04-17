import numpy as np
import matplotlib.pyplot as plt
import math
import urllib.request
import os
import time
def download(filename):
    url = "http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/"+filename
    if not os.path.exists(filename):
        urllib.request.urlretrieve(url,filename)
    else:
        print("The file %s has downloaded!"%filename)
def getdata(filename):
    data_list = []
    with open(filename) as f:
        flag = False
        while True:
            line = f.readline()
            if "EOF" in line:
                break
            elif "NODE_COORD_SECTION" in line:
                flag = True
            elif flag:
                str_list = line.strip().split(" ")
                data_list.append([float(str_list[1]),float(str_list[2])])
            else:
                continue
    return np.array(data_list)
class SA:
    def __init__(self,t0,tf,alpha,markov_len=10000):
        self.t0 = t0
        self.tf = tf
        self.alpha = alpha
        self.markov_len = markov_len
        self.sol_best = None
        self.E_best = None
        self.citys_mat = None
        self.length_list = None
        self.name = time.strftime("%Y%m%d%H%M",time.localtime(time.time()))
    def solve(self,citys_mat):
        # 获取城市的个数
        self.citys_mat = citys_mat
        length = citys_mat.shape[0]
        # 获取邻接矩阵
        citys_x = citys_mat[:,0].reshape(length,1).dot(np.ones((1,length)))
        citys_y = citys_mat[:,1].reshape(length,1).dot(np.ones((1,length)))
        citys_distance = np.sqrt(np.square(citys_x - citys_x.T) + np.square(citys_y - citys_y.T))
        E_current = np.inf #表示当前解对应的距离
        E_best = np.inf  #表示最优解
        sol_new = np.arange(0,length) #初始化解
        sol_current = sol_new.copy()
        self.sol_best = sol_new.copy()
        self.length_list = []
        t = self.t0
        epoches = 0
        while t>=self.tf:
            for step in range(self.markov_len):
                # 产生随机扰动，并产生新的解
                if np.random.random() < 0.5:
                    # 两交换
                    ind1 ,ind2 = 0,0
                    while (ind1 == ind2):
                        ind1 = math.ceil(np.random.random() * (length-1))
                        ind2 = math.ceil(np.random.random() * (length-1))
                    sol_new[ind1], sol_new[ind2] = sol_new[ind2], sol_new[ind1]
                else:
                    # 三交换
                    ind1,ind2,ind3 = 0,0,0
                    while ((ind1 == ind2) or (ind1 == ind3) or (ind2 == ind3)):
                        ind1 = math.ceil(np.random.random() * (length-1))
                        ind2 = math.ceil(np.random.random() * (length-1))
                        ind3 = math.ceil(np.random.random() * (length-1))
                    tmp_list = [ind1, ind2, ind3]
                    tmp_list = sorted(tmp_list)
                    ind1, ind2, ind3 = tmp_list

                    tmp_list1 = sol_new[ind1:ind2].copy()
                    tmp_list2 = sol_new[ind2:ind3].copy()
                    sol_new[ind1:ind1 + ind3 - ind2] = tmp_list2
                    sol_new[ind1 + ind3 - ind2:ind3] = tmp_list1
                # 计算目标函数值
                E_new = 0
                for k in range(length - 1):
                    E_new = E_new + citys_distance[sol_new[k], sol_new[k + 1]]
                E_new = E_new + citys_distance[sol_new[-1], sol_new[0]]
                if E_new < E_current:
                    E_current = E_new
                    sol_current = sol_new.copy()
                    if E_new < E_best:
                        E_best = E_new
                        self.sol_best = sol_new.copy()
                else:
                    if np.random.random()<np.exp(-(E_new-E_current)/t):
                        E_current = E_new
                        sol_current = sol_new.copy()
                    else:
                        sol_new = sol_current.copy()

                self.E_best = E_best
            self.length_list.append(E_best)
            epoches += 1
            print("epoches: %d,temperature:%.2f"%(epoches,t))
            t = t*self.alpha
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
        plt.savefig(os.path.join("log","%s_loss.png"%self.name))
        plt.close()
        # draw dots
        for k in range(0, len(self.sol_best)-1):
            start = self.citys_mat[self.sol_best[k]]
            end = self.citys_mat[self.sol_best[k+1]]
            plt.plot(start[0],start[1],"bo")
            plt.plot(end[0],end[1],"bo")
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
        plt.title(label="length:%.2f"%self.E_best)
        plt.savefig(os.path.join("log","%s_route.png"%self.name))
        plt.show()
def main():
    filename = "eil51.tsp"
    download(filename)
    data_list = getdata(filename)
    t0, tf, alpha, markov_len = 100, 3, 0.99, 10000
    model = SA(t0, tf, alpha, markov_len)
    model.solve(data_list)
    model.draw()
if __name__ == '__main__':
    main()

