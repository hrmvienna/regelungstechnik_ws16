function x = int_euler_1_2(A, B, u, x0, Ta, Te)

    % A = Dynamikmatrix
    % B = Eingangsmatrix
    % u = Eingangsgroesse
    % x0 = Anfangsbedingung
    % Ta = Abtastzeit
    % Te = Endzeitpunkt
    trange = 0:Ta:Te;     

    %Berechnung des DGL-Sys mit dem Euler-Vorwaertsverfahren:
    %Als Anfangsbedingung wird x(0) genommen und ausgehend von t=0 zur
    %gewuenschten Position (in unserem Fall length(trange)) mit der 
    %Schrittweite Ta iteriert.
    x=x0(:,1);             % : means, that all elements in the row are selected
    for k=1:length(trange)-1      %length...returns the number of elements...
        x(:,k+1) = x(:,k)+Ta*(A*x(:,k)+B*u(k)); 
    end
end