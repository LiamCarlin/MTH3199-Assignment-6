function dVdt = string_rate_func01(t,V,string_params)

n = string_params.n;
M = string_params.M;
Tf = string_params.Tf;
c = string_params.c;
dx = string_params.dx;

% Unpack the positions and velocities
U = V(1:n);
dUdt = V(n+1:2*n);

% Find the position and velocity of the moving endpoint
Uf = string_params.Uf_func(t);
dUfdt = string_params.dUfdt_func(t);

% Include both endpoints
U_full = [0; U; Uf];
dUdt_full = [0; dUdt; dUfdt];

% Calculate the acceleration of each mass
d2Udt2 = zeros(n,1);

for i = 1:n
    position_term = U_full(i) - 2*U_full(i+1) + U_full(i+2);
    velocity_term = dUdt_full(i) - 2*dUdt_full(i+1) + dUdt_full(i+2);

    d2Udt2(i) = (n/M)*(Tf/dx*position_term + c/dx*velocity_term);
end

% Assemble the state derivative
dVdt = [dUdt; d2Udt2];

end