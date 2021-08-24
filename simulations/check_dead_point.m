function [ is_dead ] = check_dead_point(x, y, A, bu, bl, Q, eps_abs)

is_dead = false;

Ax = A*x;
z = min(bu, max(Ax, bl));
eps_pri = eps_abs + eps_abs*norm([Ax;z],inf);

cons = min(abs(bu-Ax), abs(bl - Ax));

act = cons < eps_pri | abs(y) > eps_pri;

if (sum(act) == size(A, 1))
    disp('Weird, all constraints active?');
    return;
end

At = A';
Aact = At(:,act);

n = size(Aact, 2);
if n < 100000
    [~, Z, ra] = rowspaces(Aact'); Z = Z';
    nullspace_not_empty = ra < n;
else
    [~, spr] = spspaces(Aact', 2);
    J = spr{3};
    nullspace_not_empty = ~isempty(J);
    Z = spr{1};
    Z = Z(:,J);
end

if nullspace_not_empty
    H = Z'*(Q*Z); 
    H = (H+H')/2;
    lam_min = min(eig(H));
    if lam_min < -1e-10
        is_dead = true;
    end
end


end

