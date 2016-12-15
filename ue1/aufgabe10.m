
clear all;close all;clc;
% Parameter
A = [-1 2; -2 -6];

syms s

p1 = [1, 8, 19, 12] %s^3 + 8*s^2 + 19*s + 12
p2 = poly(A)
p3 = conv(p1, p2)


polydiff(p2)
