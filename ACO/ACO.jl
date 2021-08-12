using Base.Filesystem:ispath,isfile,joinpath,mkdir
using Base:download
using Printf:@sprintf
using LinearAlgebra
using Combinatorics
using Random
using StatsBase:sample
function download_dataset(root_path,filename)
    url = "http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/"*filename
    if !ispath(root_path)
        mkdir(root_path)
    end
    save_file_path = joinpath(root_path,filename)
    if !isfile(filename)
        download(url,save_file_path)
        message=@sprintf("The data set: %s downloaded!",save_file_path)
        @info message
    else
        message = @sprintf("The data set: %s already has downloaded!",save_file_path)
        @info message
    end
end
function get_data(filename)
    data_list = []
    open(filename,"r") do stream
        flag = false
        for line in eachline(stream)
            if occursin("EOF",line)
                break
            elseif occursin("NODE_COORD_SECTION",line)
                flag = true
            elseif flag
                tp_list = []
                output = split(line,['\t',' '])
                vala = parse(Int64,output[1])
                valb = parse(Float32,output[2])
                valc = parse(Float32,output[3])
                append!(tp_list,vala)
                append!(tp_list,valb)
                append!(tp_list,valc)
                append!(data_list,tp_list)
            end
        end
    end
    arr_len = Int64(length(data_list)/3)
    data_list = reshape(data_list,3,arr_len)
    return data_list'
end
function random_array(citys_num)
    tp_mat = [i for i in range(1,stop=citys_num,step=1)]
    tp_mat = shuffle(tp_mat)
    return tp_mat
end
function randseed(ant_num,citys_num)
    if ant_num <citys_num
        tp_mat = random_array(citys_num)
        initial_route = tp_mat[:ant_num]
    else
        initial_route = zeros(Int64,ant_num)
        initial_route[1:citys_num] = random_array(citys_num)
        tmp_index = citys_num+1
        tmp_left = ant_num % citys_num
        while tmp_index + citys_num <= ant_num
            initial_route[tmp_index:citys_num + tmp_index-1] = random_array(citys_num)
            tmp_index = tmp_index + citys_num
        end
        if tmp_left != 0
            tp_mat = random_array(citys_num)
            initial_route[tmp_index:end] = tp_mat[1:tmp_left]
        end
    end
    return initial_route
end

function main()
    filename = "eil51.tsp"
    root_path = "data"
    download_dataset(root_path,filename)
    data_file_path = joinpath(root_path,filename)
    citys_mat = get_data(data_file_path)
    ant_num,alpha,beta,rho,Q,epoches = 500,1,5,0.2,10,20
    # 获取城市的个数
    citys_num = size(citys_mat)[1]
    # 获取邻接矩阵
    citys_x = reshape(citys_mat[:,2],citys_num,1).*ones(1,citys_num)
    citys_y = reshape(citys_mat[:,3],citys_num,1).*ones(1,citys_num)
    citys_index = citys_mat[:,1]
    citys_distance = sqrt.((citys_x-citys_x').^2+(citys_y-citys_y').^2)
    # 初始化启发函数

    tpMatrix = zeros(citys_num,citys_num)
    for k in range(1,stop=citys_num,step=1)
        tpMatrix[k,k] = 1.0*Inf32
    end
    Heu_f = 1.0./(citys_distance + tpMatrix)
    # 信息素矩阵
    Tau_table = ones(Float32,citys_num,citys_num)
    # 每一次迭代过程中每个蚂蚁的路径记录表
    Route_table = zeros(Int64,ant_num,citys_num,)
    #Route_table = parse(Int64,Route_table)
    # 每一次迭代过程中的最佳路径
    Route_best = zeros(Int64,epoches,citys_num)
    #Route_best = parse(Int64,Route_best)
    # 每一次迭代过程中最佳路径记录表
    Length_best = zeros(Float32,epoches)
    # 每次迭代过程中蚂蚁的平均路径长度
    Length_average = zeros(Float32,epoches)
    # 每次迭代过程中当前路径长度
    Length_current = zeros(Float32,ant_num)
    iter = 1
    while iter<=epoches
        # 产生城市集合表
        # 随机产生各个蚂蚁的起点城市
        #println(size(Route_table))
        Route_table[:,1]= randseed(ant_num,citys_num)

        # 更新信息素
        Delta_tau = zeros(Float32,citys_num, citys_num)
        for k in range(1,stop=ant_num,step=1)
            # 用于记录蚂蚁下一个访问的城市集合
            # 蚂蚁已经访问过的城市
            tabu = [Route_table[k,1]]

            city_index = Route_table[k,1]
            for i=2:citys_num
                allow_set = setdiff(Set([i for i in range(1,stop=citys_num,step=1)]),Set(tabu))
                allow_set = [item for item in allow_set]
                # 初始化城市之间的转移概率
                P_table = zeros(Float32,length(allow_set))
                # 计算城市之间的转移概率

                for item in enumerate(allow_set)
                    j = item[1]
                    val = item[2]
                    P_table[j] = Tau_table[city_index,val]^alpha*Heu_f[city_index,val]^beta
                end
                P_table = P_table./sum(P_table)
                # 轮盘赌算法来选择下一个访问的城市
                #out_prob = np.cumsum(P_table)
                while true
                    prob_r = rand()
                    index_need = findall(x->x >= prob_r,P_table)
                    if length(index_need)>0
                        city_index2 = allow_set[index_need[1]]
                        Route_table[k,i] = city_index2
                        append!(tabu,city_index2)
                        city_index = city_index2
                        break
                    end
                end
            end
            append!(tabu,tabu[1])
            # 计算蚂蚁路径的距离信息
            for j=1:citys_num
                Length_current[k] = Length_current[k] + citys_distance[tabu[j],tabu[j+1]]
            end
            for j in range(1,stop=citys_num,step=1)
                Delta_tau[tabu[j],tabu[j+1]] = Delta_tau[tabu[j],tabu[j+1]] + Q / Length_current[k]
            end
        end
        #println(typeof(Route_table))
        #println(size(Route_table))
        # 计算最短路径、最短路径长度以及平均路径长度
        minval = findmin(Length_current)
        Length_best[iter],index = minval[1],minval[2]
        Route_best[iter,:] = Route_table[index,:]

        Length_average[iter] = sum(Length_current)/length(Length_current)
        #更新信息素
        Tau_table = (1-rho).*Tau_table + Delta_tau
        #Route_table = np.zeros((self.ant_num,citys_num),dtype=np.int)
        Length_current = zeros(ant_num)
        message = @sprintf("epoches:%d,best value every epoches is %.4f",iter,Length_best[iter])
        @info message
        iter += 1
        Array
    end
end
main()
