function [B0,L,A,E] = elementInfo(e, IX, X, mprop)

    % area and youngs
    propno = IX(e,3);
    A = mprop(propno,2);
    E = mprop(propno,1);
    
    % dX,dY and length
    nodeNumber1 = IX(e,1);
    nodeNumber2 = IX(e,2);
    
    node1 = X(nodeNumber1,:);
    node2 = X(nodeNumber2,:);
    
    dX = node2(1)-node1(1);
    dY = node2(2)-node1(2);
    
    L = sqrt(dX^2+dY^2); % used later to calc total length
    
    % B0 (original geometry) 
    B0 = (1/(L)^2) * [-dX ; -dY;  dX;  dY];
    