% Basis truss program
% Modified 3/9-2006 by OS
% Modified 20/8-2012 by JSJ
clear all
close all
clear
clf
example41;
%example1;

tic;
timerVal = tic;

neqn = size(X,1)*size(X,2);         % Number of equations
ne = size(IX,1);                    % Number of elements
disp(['Number of DOF ' sprintf('%d',neqn) ...
    ' Number of elements ' sprintf('%d',ne)]);

% Initialize Global Variables
K=zeros(neqn,neqn);                     % Stiffness matrix
P=zeros(neqn,1);                        % Force vector
D=zeros(neqn,1);                        % Displacement vector
R=zeros(neqn,1);                        % Residual vector
strain=zeros(ne,1);                     % Element strain vector
stress=zeros(ne,1);                     % Element stress vector

V = 6;                                  % maximum volume
v = zeros(ne,1);                        % element volumes
p = 1.5;                                % power of penalization
dg = zeros(ne,1);
df = zeros(ne,1);
epsilon = 10^(-3);                      % breaking condition
max_iopt = 1000;
rho = ones(ne,1);
rho_min = 10^(-6);
sed = zeros(ne,1);
ss = zeros(ne,1);

% calculate element volumes
for e = 1: ne
    [B0, L, A, E] = elementInfo(e,IX,X,mprop);
    v(e) = L*A;
end

    % the dg is given by:
    dg = v;

% sets rho to satisfi g = 0 -> v'rho = V
value_rho = V/sum(v);
rho = rho*value_rho;

%%
for iopt = 1:max_iopt
    
   %% 
    rho_old = rho;
        
    %% K MATRIX
    K=zeros(neqn,neqn);                     % Stiffness matrix
    
    for e = 1 : ne
        
        % pull information from element
        [B0, L, A, E] = elementInfo(e,IX,X, mprop);
        
        % create the element stiffness matrix
        ke = createElementK( A, E, L, B0);
        
        % create edof vector of element
        edof = createEdof(e, IX);
        
        % assembly of the total stiffness matrix
        K( edof, edof ) = K(edof, edof) + rho(e)*ke;
        
        
    end
    
    %% BOUNDARY CONDITIONS
    
    % Applying loads to P vector
    for ind1 = 1:size(loads,1)
        
        P(2*loads(ind1,1)- (2-loads(ind1,2))) = loads(ind1,3);
        
    end
    
    % Apply boundary conditions and removing stuff
    for ind2 = 1:size(bound,1)
        
        % for simplicity defining a vector bc with indexes
        bc(ind2) = 2*bound(ind2,1)-(2-bound(ind2,2));
        
        % sets rows and cols = 0
        K(:,bc(ind2)) = 0; %
        K(bc(ind2) ,:) = 0;
        % ... and the diagonal element = 1
        K(bc(ind2),bc(ind2)) = 1;
        
        % sets the load to equal the displacement
        P(bc(ind2)) = 0 ;
        
    end
    
    %%
    % calculating the global displacement vector
    D = K\P;
    
    % total length of structure
    totalLength = sum(L);
    
    % displacement at loads
    dispLoad = [ D(2*loads(:,1)) D(2*loads(:,1)-1) ];
     
    f(iopt) = D'*P; 

%%   CALCULATING df AND dg

for e = 1:ne
    
    % pull information from element
    [B0, L, A, E] = elementInfo(e,IX,X, mprop);
    
    % create the element stiffness matrix
    ke = createElementK( A, E, L, B0);
    
    % create edof vector of element
    edof = createEdof(e, IX);
    
    % disp
    d = D(edof);
    
    % the df is given by:
    df(e) = -p*rho(e)^(p-1)*d'*ke*d;
        
end

%% CALCULATING THE RHO
[ rho, lambda] = bisect( rho_old, rho_min, V, df, dg, ne);

 
%% BREAKING CONDITION
if norm(rho_old - rho) < epsilon * norm(rho)
    break;
end

end

%% CALCULATING STRESS STRAIN DENSITY
% calculating stress and strain of the elements
for e=1:ne
    
    % pull information from element
    [B0, L, A, E] = elementInfo(e,IX,X, mprop);
    
    % create edof vector of element
    edof = createEdof(e, IX);
    
    % create the displacement vector of the element
    d = D(edof);
    
    % calculate strain and stress
    strain(e) = B0'*d;
    stress(e) = rho(e)^p*E*B0'*d;
    ss(e,1) = strain(e)*stress(e);
    
    % checking the strain energy density of elements not on bounds
    if rho(e) > rho_min
        if rho(e) < 1
        sed(e) = ss(e)/rho(e);
        end
    end
   
    % 
figure(2) 
clf
plot(1:iopt, f) 
legend('Compliance vs itterations')
xlabel('Itterations of topology optimization') % x-axis label
ylabel('Compliance') % y-axis label

% Plotting Un-Deformed and Deformed Structure
figure(3)
clf
hold on

for e = 1:ne
    xx = X(IX(e,1:2),1);
    yy = X(IX(e,1:2),2);
    plot(xx,yy,'k:','LineWidth',1.)
    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];
    xx = xx + D(edof(1:2:4));
    yy = yy + D(edof(2:2:4));
    
    if  rho(e) > 10*rho_min
        plot(xx,yy,'b','LineWidth',5*rho(e))
    end
    
end
axis equal;
plotsupports
plotloads
hold off
legend('Truss structure')

    
    
end

toc
fprintf('Topology itterations stopped at iopt = %6f \n',iopt)


figure(1)
bar(sed)
legend('Strain energy density of elements not on bounds')
xlabel('Elements') % x-axis label
ylabel('Energy') % y-axis label


% 
figure(2) 
plot(1:iopt, f) 
legend('Compliance vs itterations')
xlabel('Itterations of topology optimization') % x-axis label
ylabel('Compliance') % y-axis label

% Plotting Un-Deformed and Deformed Structure
figure(3)
hold on

for e = 1:ne
    xx = X(IX(e,1:2),1);
    yy = X(IX(e,1:2),2);
    plot(xx,yy,'k:','LineWidth',1.)
    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];
    xx = xx + D(edof(1:2:4));
    yy = yy + D(edof(2:2:4));
    
    if  rho(e) > 10*rho_min
        plot(xx,yy,'b','LineWidth',5*rho(e))
    end
    
end
axis equal;
plotsupports
plotloads
hold off
legend('Truss structure')
