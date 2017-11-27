function [ ke ] = createElementK( A, E, L, B0)
    % Element stiffness
    ke = A*E*L*B0*B0';
end

