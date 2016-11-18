% ===========
% Aufgabe 1.1:
% ===========

clear all   % Es werden alle Variablen gelöscht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zurückgesetzt

%% 
clear all;close all;clc;

% Aufgabe 1.1.1:

% Matrix A und Vektor b
A = [1, 2, 4; 2, 2, 1; 3, 2, 0];
b = [5; 4; 1];

% Zeit
tic
x1 = inv(A)*b
toc
tic
x2 = mldivide(A, b)
toc

% numerischer Fehler
norm(A*x1 - b)
norm(A*x2 - b)

% mldivide ist schneller (0.000148 vs 0.000070) und numerisch genauer

%% 
clear all;close all;clc;

% Aufgabe 1.1.2:

% rect(x) ~~ A* sum(k=1 bis n) 4/pi*(2k-1) * sin((2k-1)x)

A = 10;
x = 0:1:10;
n = 1;
y = 0

for n=1:1:100
    y = y + A*(4/(pi*(2*n -1))) * sin((2*n -1)*x);
    plot(y)
    title(['rec(x) for n = ', num2str(n)])
    grid on
    pause()
    %pause(0.01)
end

%% 
clear all;close all;clc;

% Aufgabe 1.1.3:

% f(x,y,t) = sin(x*pi/10)*sin(y*pi/10)*|sin(t)|
step_size = 0.2;
range = 0:step_size:10;
times = [0.1, 0.5, 1, 2, 10, 32];

[xx,yy] = meshgrid(range,range); % Koordinatenmatrix erzeugen
figure
for t = times

    zz = sin((xx*pi)/10).*sin((yy*pi)/10).*abs(sin(t));

    % Gitterplot
    mesh(xx,yy,zz)
    box on
    grid on
    
    suptitle(['t = ', num2str(t)]);
    pause()
end

%%
% ===========
% Aufgabe 1.2:
% ===========

%% 
clear all;close all;clc;

% Aufgabe 1.2.1
