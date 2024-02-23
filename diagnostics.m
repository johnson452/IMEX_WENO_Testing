function [app] = diagnostics(app)
% Creates the diagnostics for our problem

clf()

% Create plots:
subplot(2,3,1)
f = app.f_IC;
[X,V] = meshgrid(app.grid.x,app.grid.v);
contourf(X',V',f)
title("f(x,v) (I.C.)")
xlabel("x")
ylabel("v")


subplot(2,3,2)
f = app.f;
[X,V] = meshgrid(app.grid.x,app.grid.v);
contourf(X',V',f)
title("f(x,v)")
xlabel("x")
ylabel("v")

subplot(2,3,3)
plot(app.grid.x,app.n)
hold on
plot(app.grid.x,app.n_midpoint1)
hold on
plot(app.grid.x,app.n_Eq)
hold on
plot(app.grid.x,app.n_IC)
title("n(x)")
xlabel("x")
ylabel("n")
legend("n(t)","n_{mid}","n_{Eq}","n_{I.C.}")

subplot(2,3,4)
plot(app.grid.x,app.u)
hold on
plot(app.grid.x,app.u_midpoint1)
hold on
plot(app.grid.x,app.u_Eq)
hold on
plot(app.grid.x,app.u_IC)
title("u(x)")
xlabel("x")
ylabel("u")
legend("u(t)","u_{mid}","u_{Eq}","u_{I.C.}")

subplot(2,3,5)
plot(app.grid.x,app.T)
hold on
plot(app.grid.x,app.T_midpoint1)
hold on
plot(app.grid.x,app.T_Eq)
hold on
plot(app.grid.x,app.T_IC)
title("T(x)")
xlabel("x")
ylabel("T")
legend("T(t)","T_{mid}","T_{Eq}","T_{I.C.}")


pause(0.1)

end