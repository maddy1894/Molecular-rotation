f = @(t,y)[2*y(1)+3*y(2) ; 3*y(1)+2*y(2)]; % define differential equations
tspan = [0 20];                             % set interval of integration
y0=[0,1]
[T Y] = ode45(f,tspan,y0);                                    % call integrator
plot(T,Y(:,1),'-o',T,Y(:,2),'-o')
display(Y)
