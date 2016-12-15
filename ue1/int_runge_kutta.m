function x = int_runge_kutta(A, B, u, x0, Ta, Te)

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
        k1 = (A*x(:,k)+B*u(k));
        
        uk_p_halb = (u(k) + u(k+1))/2;
        
        k2 = (A*(x(:,k) + (1/2)*Ta*k1)+B*uk_p_halb);
        
        k3 = (A*(x(:,k) + (1/2)*Ta*k2) + B*uk_p_halb);
        
        k4 = (A*(x(:,k) + Ta * k3) + B*(u(k+1)));
        
        x(:,k+1) = x(:,k) + (1/6)*Ta*(k1 + 2*k2 + 2*k3 + k4); 
    end
end