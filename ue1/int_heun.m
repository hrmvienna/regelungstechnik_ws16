function x = int_heun(A, B, u, x0, Ta, Te)

    % A = Dynamikmatrix
    % B = Eingangsmatrix
    % u = Eingangsgroesse
    % x0 = Anfangsbedingung
    % Ta = Abtastzeit
    % Te = Endzeitpunkt
    trange = 0:Ta:Te;     
    
    s0 = 1;
    vo= 0.5;
    m = 1;
    c = 1;
    k = 2;

    %Berechnung des DGL-Sys mit dem Euler-Vorwaertsverfahren:
    %Als Anfangsbedingung wird x(0) genommen und ausgehend von t=0 zur
    %gewuenschten Position (in unserem Fall length(trange)) mit der 
    %Schrittweite Ta iteriert.
    x=x0(:,1);             % : means, that all elements in the row are selected
    
    xkp = x0(:,1);
    
    for k=1:length(trange)-1      %length...returns the number of elements...
        xkp(:,k+1) = x(:,k)+Ta*(A*x(:,k)+B*u(k));
        
        x(:,k+1) = (1/2)* x(:,k)+(1/2)*(xkp(:,k) + Ta*(A*xkp(:,k)+B*u(k+1))); 
    end
end