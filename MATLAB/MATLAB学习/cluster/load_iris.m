clc;clear;
[A,B,C,D,E] = textread('iris.data','%f%f%f%f%s','delimiter',',');
Data=[A,B,C,D];
Group=D;