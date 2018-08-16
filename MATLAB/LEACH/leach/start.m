clc, clear all, close all

numNodes = 100; % number of nodes
p = 0.1;

netArch  = newNetwork(100, 100, 50, 175); %%%%%%%%%%   网络范围     %%%%%%%%%
                                          % NetArch.Yard.Type           网络类型：矩形
                                          % NetArch.Yard.Length         网络长：100
                                          % NetArch.Yard.Width          网络宽：100
                                          % NetArch.Sink.x              sink（X）：50
                                          % NetArch.Sink.y              sink（Y）：175
                                         %%%%%%%%%%%  能量模型     %%%%%%%%
                                          % NetArch.Energy.init         0.5
                                          % NetArch.Energy.transfer     50*0.000000001
                                          % NetArch.Energy.receive      50*0.000000001
                                          % NetArch.Energy.freeSpace    10*0.000000000001
                                          % NetArch.Energy.multiPath    0.0013*0.000000000001
                                          % NetArch.Energy.aggr         5*0.000000001
                                          
nodeArch = newNodes(netArch, numNodes);   % nodeArch.nodesLoc(i, 1)，nodeArch.nodesLoc(i, ) 设置节点参数——GPS坐标
                                          %%%%%%%%%%%%    节点模型     %%%%%%%
                                          % nodeArch.node(i).G          0
                                          % nodeArch.node(i).type       'N'
                                          % nodeArch.node(i).energy     0.5
                                          % nodeArch.node(i).CH         -1
                                          % nodeArch.dead(i)            0
                                          %%%%%%%%%%%%   生存、死亡节点数量 %%%%%
                                          % nodeArch.numNode  
                                          % nodeArch.numDead 
                                          
roundArch = newRound();% NetRound.numRound     最大轮数 9999
                       % NetRound.packetLength   CH 到 BH 数据报  6400 
                       % NetRound.ctrPacketLength    node 到 CH 数据报 200 

plot1

par = struct;

for r = 1:roundArch.numRound
    r
    clusterModel = newCluster(netArch, nodeArch, 'leach', r, p);
    clusterModel = dissEnergyCH(clusterModel, roundArch);
    clusterModel = dissEnergyNonCH(clusterModel, roundArch);
    nodeArch     = clusterModel.nodeArch; % new node architecture after select CHs
    
    par = plotResults(clusterModel, r, par);
    if nodeArch.numDead == nodeArch.numNode
        break
    end
end


