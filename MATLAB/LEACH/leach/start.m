clc, clear all, close all

numNodes = 100; % number of nodes
p = 0.1;

netArch  = newNetwork(100, 100, 50, 175); %%%%%%%%%%   ���緶Χ     %%%%%%%%%
                                          % NetArch.Yard.Type           �������ͣ�����
                                          % NetArch.Yard.Length         ���糤��100
                                          % NetArch.Yard.Width          �����100
                                          % NetArch.Sink.x              sink��X����50
                                          % NetArch.Sink.y              sink��Y����175
                                         %%%%%%%%%%%  ����ģ��     %%%%%%%%
                                          % NetArch.Energy.init         0.5
                                          % NetArch.Energy.transfer     50*0.000000001
                                          % NetArch.Energy.receive      50*0.000000001
                                          % NetArch.Energy.freeSpace    10*0.000000000001
                                          % NetArch.Energy.multiPath    0.0013*0.000000000001
                                          % NetArch.Energy.aggr         5*0.000000001
                                          
nodeArch = newNodes(netArch, numNodes);   % nodeArch.nodesLoc(i, 1)��nodeArch.nodesLoc(i, ) ���ýڵ��������GPS����
                                          %%%%%%%%%%%%    �ڵ�ģ��     %%%%%%%
                                          % nodeArch.node(i).G          0
                                          % nodeArch.node(i).type       'N'
                                          % nodeArch.node(i).energy     0.5
                                          % nodeArch.node(i).CH         -1
                                          % nodeArch.dead(i)            0
                                          %%%%%%%%%%%%   ���桢�����ڵ����� %%%%%
                                          % nodeArch.numNode  
                                          % nodeArch.numDead 
                                          
roundArch = newRound();% NetRound.numRound     ������� 9999
                       % NetRound.packetLength   CH �� BH ���ݱ�  6400 
                       % NetRound.ctrPacketLength    node �� CH ���ݱ� 200 

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


