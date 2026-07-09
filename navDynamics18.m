%NAVIGATION/STATE ESTIMATION MODEL

%STATES: xE, yE, zE, u, v, w, phi, theta, psi, lambda_Ax, lambda_Ay, lambda_Az, lambda_p, lambda_q, lambda_r, Wx, Wy, Wz

function xdot = navDynamics18(x, u_m)
    g = 9.80665;

    % States
    xE    = x(1);
    yE    = x(2);
    zE    = x(3);

    u     = x(4);
    v     = x(5);
    w     = x(6);

    phi   = x(7);
    theta = x(8);
    psi   = x(9);

    lambda_Ax   = x(10);
    lambda_Ay   = x(11);
    lambda_Az   = x(12);

    lambda_p    = x(13);
    lambda_q    = x(14);
    lambda_r    = x(15);

    Wx    = x(16);
    Wy    = x(17);
    Wz    = x(18);

    % Measured IMU inputs
    Ax_m = u_m(1);
    Ay_m = u_m(2);
    Az_m = u_m(3);

    p_m  = u_m(4);
    q_m  = u_m(5);
    r_m  = u_m(6);

    % Bias-corrected IMU inputs
    Ax = Ax_m - lambda_Ax;
    Ay = Ay_m - lambda_Ay;
    Az = Az_m - lambda_Az;

    p = p_m - lambda_p;
    q = q_m - lambda_q;
    r = r_m - lambda_r;

    % Useful trigonometric terms
    cphi = cos(phi);
    sphi = sin(phi);

    ctheta = cos(theta);
    stheta = sin(theta);
    ttheta = tan(theta);

    cpsi = cos(psi);
    spsi = sin(psi);

    % State Derivatives
    xE_dot = (u*ctheta + (v*sphi + w*cphi)*stheta)*cpsi ...
             - (v*cphi - w*sphi)*spsi ...
             + Wx;

    yE_dot = (u*ctheta + (v*sphi + w*cphi)*stheta)*spsi ...
             + (v*cphi - w*sphi)*cpsi ...
             + Wy;

    zE_dot = -u*stheta ...
             + (v*sphi + w*cphi)*ctheta ...
             + Wz;

    % Body velocity dynamics
    u_dot = Ax - g*stheta + r*v - q*w;
    v_dot = Ay + g*sphi*ctheta + p*w - r*u;
    w_dot = Az + g*cphi*ctheta + q*u - p*v;

    % Attitude dynamics
    phi_dot   = p + q*sphi*ttheta + r*cphi*ttheta;
    theta_dot = q*cphi - r*sphi;
    psi_dot   = (q*sphi + r*cphi)/ctheta;

    % Biases and wind are modeled as constants
    bAx_dot = 0;
    bAy_dot = 0;
    bAz_dot = 0;

    bp_dot = 0;
    bq_dot = 0;
    br_dot = 0;

    Wx_dot = 0;
    Wy_dot = 0;
    Wz_dot = 0;

    % Full state derivative
    xdot = [xE_dot;
            yE_dot;
            zE_dot;
            u_dot;
            v_dot;
            w_dot;
            phi_dot;
            theta_dot;
            psi_dot;
            bAx_dot;
            bAy_dot;
            bAz_dot;
            bp_dot;
            bq_dot;
            br_dot;
            Wx_dot;
            Wy_dot;
            Wz_dot];
end