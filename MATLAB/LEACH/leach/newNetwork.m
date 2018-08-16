function NetArch = newNetwork(Length, Width, sinkX, sinkY, initEnergy...
    , transEnergy, recEnergy, fsEnergy, mpEnergy, aggrEnergy)
% Create the network architecture with desired parameters
%   
%   Input:
%       Length      Length of the yard
%       Width       Width of the yard
%       sinkX       x cordination of base station
%       sinkY       y cordination of base station
%       initEnergy  Initial energy of each node
%       transEnergy Enery for transferring of each bit (ETX)
%       recEnergy   Enery for receiving of each bit (ETX)
%       fsEnergy    Energy of free space model
%       mpEnergy    Energy of multi path model
%       aggrEnergy  Data aggregation energy     
%   Example:
%       NetArch = createNetwork();
%
% Hossein Dehghan, hd.dehghan@gmail.com
% Ver 1. 2/2013

    %%%% Create the yard
    NetArch.Yard.Type = 'Rect'; % Rectangular
    if ~exist('Length','var')
        NetArch.Yard.Length = 100; % default of the yard is 100 in x cordination
    else
        NetArch.Yard.Length = Length;
    end
    if ~exist('Width','var')
        NetArch.Yard.Width = 100; % default of the yard is 100 in y cordination
    else
        NetArch.Yard.Width = Width;
    end
    
    %%%% Create base station
    % x and y Coordinates of the base station
    % default of the base station is in the center of the yard
    if ~exist('sinkX','var')
        NetArch.Sink.x = Yard.Length / 2;
    else
        NetArch.Sink.x = sinkX;
    end
    if ~exist('sinkY','var')
        NetArch.Sink.y = Yard.Width / 2;
    else
        NetArch.Sink.y = sinkY;
    end

    %%%% Energy Model (all values in Joules)
    % Initial Energy
    if ~exist('initEnergy','var')
        NetArch.Energy.init = 0.5; 
    else
        NetArch.Energy.init = initEnergy; 
    end
    
    % Enery for transferring of each bit (ETX)
    if ~exist('transEnergy','var')
        NetArch.Energy.transfer = 50*0.000000001;
    else
        NetArch.Energy.transfer = transEnergy; 
    end
    if ~exist('recEnergy','var')
        NetArch.Energy.receive = 50*0.000000001;
    else
        NetArch.Energy.receive = recEnergy; 
    end
    
    % Transmit Amplifier types
    if ~exist('recEnergy','var')
        NetArch.Energy.freeSpace = 10*0.000000000001;
    else
       NetArch.Energy.freeSpace = fsEnergy; 
    end
    if ~exist('recEnergy','var')
        NetArch.Energy.multiPath = 0.0013*0.000000000001;
    else
        NetArch.Energy.multiPath = mpEnergy; 
    end
    
    %Data Aggregation Energy
    if ~exist('recEnergy','var')
        NetArch.Energy.aggr = 5*0.000000001;
    else
        NetArch.Energy.aggr = aggrEnergy; 
    end


end