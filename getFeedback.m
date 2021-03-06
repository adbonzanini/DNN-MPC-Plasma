function [u_mpc, feas, v_opt, args] = getFeedback(solver, args, x, ytarget)

% Extract arguments
w0 = args.w0;
lbw = args.lbw;
ubw = args.ubw;
lbg = args.lbg;
ubg = args.ubg;
offset_mpc = args.offset_mpc;
offset_x0 = args.offset_x0;
offset_ytar = args.offset_ytar;
nx = args.nx;
ny = args.ny;
nu = args.nu;

% Update initial state
w0(offset_x0+1:offset_x0+nx) = x;
lbw(offset_x0+1:offset_x0+nx) = x;
ubw(offset_x0+1:offset_x0+nx) = x;

% Update target
w0(offset_ytar+1:offset_ytar+ny) = ytarget;
lbw(offset_ytar+1:offset_ytar+ny) = ytarget;
ubw(offset_ytar+1:offset_ytar+ny) = ytarget;

% Solve the nmpc problem
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);
u_mpc = w_opt(offset_mpc+1:offset_mpc+nu);

% Check feasibility
ipopt_stats = solver.stats();
if strcmp(ipopt_stats.return_status,'Solve_Succeeded')
    feas = 1;
else
    feas = -1;
end

% Get objective value
v_opt = full(sol.f);

% Warm start, if specified
if args.warm_start == 1
    args.w0 = w_opt;
end

end