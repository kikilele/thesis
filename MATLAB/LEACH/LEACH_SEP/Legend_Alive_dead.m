 rounds=1:1:301;
close all
% plotting the graph for no: of dead nodes

figure
plot(rounds,DEAD_leach,'r',rounds,DEAD_SEP)
hleg1 = legend('Leach','SEP');
xlabel('Number of Rounds')
ylabel('No: of Nodes (DEAD)')
grid on
axis tight

% plotting the graph for no: of alive nodes
figure
plot(rounds,ALIVE_Leach,'r',rounds,Live_SEP)
hleg1 = legend('Leach','SEP');
xlabel('Number of Rounds')
ylabel('No: of Nodes (ALIVE)')
grid on
axis tight