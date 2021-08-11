using Base.Filesystem:ispath
using Printf:@sprintf
using Plots
function download_dataset(filename)
    url = "http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/"*filename
    #println(url)
    if !isfile(filename)
        download(url,filename)
    else
        message = @sprintf("The file %s is exists!",filename)
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
                tp_list = Int64[]
                output = split(line,['\t',' '])
                vala = parse(Float32,output[1])
                valb = parse(Float32,output[2])
                valc = parse(Float32,output[3])
                append!(tp_list,vala)
                append!(tp_list,valb)
                append!(tp_list,valc)
                append!(data_list,tp_list)
            end
        end
        return data_list
    end
    arr_len = Int64(length(data_list)/3)
    data_list = reshape(data_list,3,arr_len)
    return data_list'
end
function main()
    filename = "eil51.tsp"
    download_dataset(filename)
    citys_mat = get_data(filename)
    # 初始化参数
    t0, tf, alpha, markov_len = 100, 1, 0.99, 10000
    # 获取城市的个数
    citys_num = size(citys_mat)[1]
    # 获取邻接矩阵
    citys_x = reshape(citys_mat[:,2],citys_num,1).*ones(1,citys_num)
    citys_y = reshape(citys_mat[:,3],citys_num,1).*ones(1,citys_num)
    citys_index = citys_mat[:,1]
    citys_distance = sqrt.((citys_x-citys_x').^2+(citys_y-citys_y').^2)
    E_current = Inf32 # 当前解对应的距离
    E_best = Inf32 # 表示的是最优解
    sol_new = [i for i in range(1,stop=citys_num,step=1)] # 初始化解
    sol_current =copy(sol_new)
    sol_best = copy(sol_new)
    len_list = Float32[]
    t = t0
    epoches = 0
    while t>=tf
        for step=1:markov_len
            # 产生随机扰动
            if rand() < 0.5
                # 两交换
                ind1,ind2 = 0,0
                while (ind1 == ind2)
                    ind1 = Int(ceil(rand()*citys_num))
                    ind2 = Int(ceil(rand()*citys_num))
                    sol_new[ind1],sol_new[ind2] = sol_new[ind2],sol_new[ind1]
                end
            else
                ind1,ind2,ind3 = 0,0,0
                while (ind1 == ind2) || (ind1==ind3)|| (ind2==ind3)
                    ind1 = Int(ceil(rand()*citys_num))
                    ind2 = Int(ceil(rand()*citys_num))
                    ind3 = Int(ceil(rand()*citys_num))
                end
                tp_list = [ind1,ind2,ind3]
                tp_list = sort(tp_list)
                ind1,ind2,ind3 = tp_list[1],tp_list[2],tp_list[3]
                tpa = copy(sol_new[ind1:ind2-1])
                tpb = copy(sol_new[ind2:ind3-1])
                sol_new[ind1:ind1+ind3-ind2-1],sol_new[ind1+ind3-ind2:ind3-1] = tpb,tpa
            end
            # 计算目标函数值
            E_new = 0
            for k in range(1,stop=citys_num-1,step=1)
                E_new+= citys_distance[sol_new[k],sol_new[k+1]]
            end
            E_new += citys_distance[sol_new[1],sol_new[end]]
            if E_new <E_current
                E_current = E_new
                sol_current = copy(sol_new)
                if E_new < E_best
                    E_best = E_new
                    sol_best = copy(sol_new)
                end
            else
                if rand()<exp(-(E_new-E_current)/t)
                    E_current = E_new
                    sol_current = copy(sol_new)
                else
                    sol_new = copy(sol_current)
                end
            end
        end
        append!(len_list,E_best)
        epoches += 1
        message = @sprintf("epoches: %d,temperature:%.2f,length:%0.4f",epoches,t,E_best)
        if E_best<0.1
            break
        end
        t = t*alpha
        @info message
    end
    #draw(len_list,sol_best,citys_mat)
end
#=
function draw(len_list,points_list,citys_mat)
    #plot(len_list,linewidth=2,title="The best length")
    #savefig("length.png")
    #closeall()
    len_arr = length(points_list)
    for k in range(1,stop=len_arr-1,step=1)
        start_id = points_list[k]
        end_id = points_list[k+1]
        tp_x = [citys_mat[start_id,2],citys_mat[end_id,2]]
        tp_y = [citys_mat[start_id,3],citys_mat[end_id,3]]
        plot!(tp_x,tp_y,linewidth=2)
    end
    start_id = points_list[end]
    end_id = points_list[1]
    tp_x = [citys_mat[start_id,2],citys_mat[end_id,2]]
    tp_y = [citys_mat[start_id,3],citys_mat[end_id,3]]
    title!("The path of the route")
    savefig("points.png")
    closeall()
end
=#
main()
