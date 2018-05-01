function Y = RK2_integrator(odefun,tspan,y0,varargin)
% RK2 integration routine 
% AUTHOR : Haresh Karnan
% Date : 30/4/2018
% Contact : haresh.miriyala@gmail.com
if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. '];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  
RK2a = 0.5;
RK2b = 0.5;
neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);
F = zeros(neq,4);

Y(:,1) = y0;
for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1);
  yi = Y(:,i-1);
  k1 = hi*feval(odefun,ti,yi,varargin{:});
  k2 = hi*feval(odefun,ti+hi,yi+k1,varargin{:});
  Y(:,i) = yi + (RK2a*k1+RK2b*k2); % second order Runge-Kutta Method (RK2)
end

Y = Y.';