%OBSERVATION MODEL FOR THE PREDICTED OUTPUTS

%MEASURED OUTPUT: x_GPS, y_GPS, z_GPS, u_GPS, v_GPS, w_GPS, phi_GPS,
%theta_GPS, psi_GPS, V_m, alpha_m, beta_m

%PREDICTION OBSERVATION MODEL

function zhat = measurementModel12(x)

xE    = x(1);
yE    = x(2);
zE    = x(3);

u     = x(4);
v     = x(5);
w     = x(6);

phi   = x(7);
theta = x(8);
psi   = x(9);

Wx    = x(16);
Wy    = x(17);
Wz    = x(18);

% Trigonometric terms
cphi = cos(phi);
sphi = sin(phi);

ctheta = cos(theta);
stheta = sin(theta);

cpsi = cos(psi);
spsi = sin(psi);

% Predicted GPS position
x_GPS_hat = xE;
y_GPS_hat = yE;
z_GPS_hat = zE;

% Predicted GPS ground-speed components
u_GPS_hat = (u*ctheta + (v*sphi + w*cphi)*stheta)*cpsi ...
    - (v*cphi - w*sphi)*spsi ...
    + Wx;

v_GPS_hat = (u*ctheta + (v*sphi + w*cphi)*stheta)*spsi ...
    + (v*cphi - w*sphi)*cpsi ...
    + Wy;

w_GPS_hat = -u*stheta ...
    + (v*sphi + w*cphi)*ctheta ...
    + Wz;

% Predicted GPS attitude
phi_GPS_hat   = phi;
theta_GPS_hat = theta;
psi_GPS_hat   = psi;

% Predicted airdata
V_hat     = sqrt(u^2 + v^2 + w^2);
alpha_hat = atan2(w, u);
beta_hat  = atan2(v, sqrt(u^2 + w^2));

% Full 12-output measurement model
zhat = [x_GPS_hat;
    y_GPS_hat;
    z_GPS_hat;
    u_GPS_hat;
    v_GPS_hat;
    w_GPS_hat;
    phi_GPS_hat;
    theta_GPS_hat;
    psi_GPS_hat;
    V_hat;
    alpha_hat;
    beta_hat];
end