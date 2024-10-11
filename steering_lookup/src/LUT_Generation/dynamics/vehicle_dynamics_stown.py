import math
import math
import casadi as ca

__author__ = "Matthias Althoff"
__copyright__ = "TUM Cyber-Physical Systems Group"
__version__ = "2020a"
__maintainer__ = "Gerald Würsching"
__email__ = "commonroad@lists.lrz.de"
__status__ = "Released"

def vehicle_dynamics_st_delayed(x, uInit, p, type):
    """
    vehicleDynamics_st - single-track vehicle dynamics
    reference point: center of mass

    Syntax:
        f = vehicleDynamics_st(x,u,p)

    Inputs:
        :param x: vehicle state vector
        :param uInit: vehicle input vector
        :param p: vehicle parameter vector

    Outputs:
        :return f: right-hand side of differential equations

    Author: Matthias Althoff
    Written: 12-January-2017
    Last update: 16-December-2017
                 03-September-2019
    Last revision: 17-November-2020
    """

    #------------- BEGIN CODE --------------

    # set gravity constant
    g = 9.81  #[m/s^2]

    #create equivalent bicycle parameters
    mu = p.mu
    if type == "pacejka":
        B_f = p.C_Pf[0]
        C_f = p.C_Pf[1]
        D_f = p.C_Pf[2]
        E_f = p.C_Pf[3]
        B_r = p.C_Pr[0]
        C_r = p.C_Pr[1]
        D_r = p.C_Pr[2]
        E_r = p.C_Pr[3]
    elif type == "linear":
        C_Sf = p.C_Sf #-p.tire.p_ky1/p.tire.p_dy1  
        C_Sr = p.C_Sr #-p.tire.p_ky1/p.tire.p_dy1  
    lf = p.l_f
    lr = p.l_r
    h = p.h_cg 
    m = p.m 
    I = p.I_z
    tau_steer = p.tau_steer

    #states
    #x0 = x-position in a global coordinate system
    #x1 = y-position in a global coordinate system
    #x2 = steering angle of front wheels
    #x3 = yaw angle
    #x4 = velocity in x-direction
    #x5 = velocity in y direction
    #x5 = yaw rate

    #u0 = steering angle
    #u1 = longitudinal acceleration

    u = uInit

    # system dynamics

    # compute lateral tire slip angles
    alpha_f = -math.atan((x[5] + x[6] * lf) / x[4]) + x[2] 
    alpha_r = -math.atan((x[5] - x[6] * lr) / x[4])


    # compute vertical tire forces
    F_zf = m * (-u[1] * h + g * lr) / (lr + lf)
    F_zr = m * (u[1] * h + g * lf) / (lr + lf)

    F_yf = F_yr = 0

    # combined slip lateral forces
    if type == "pacejka":
        F_yf = mu * F_zf * D_f * math.sin(C_f * math.atan(B_f * alpha_f - E_f*(B_f * alpha_f - math.atan(B_f * alpha_f))))
        F_yr = mu * F_zr * D_r * math.sin(C_r * math.atan(B_r * alpha_r - E_r*(B_r * alpha_r - math.atan(B_r * alpha_r))))
    elif type == "linear":
        F_yf = mu * F_zf * C_Sf * alpha_f
        F_yr = mu * F_zr * C_Sr * alpha_r

    f = [x[4]*math.cos(x[3]) - x[5]*math.sin(x[3]), 
        x[4]*math.sin(x[3]) + x[5]*math.cos(x[3]), 
        (u[0] - x[2])/0.2,
        x[6], 
        u[1], 
        1/m * (F_yr + F_yf) - x[4] * x[6],
        1/I * (-lr * F_yr + lf * F_yf)] 
        # -mu*m/(x[4]*I*(lr+lf))*(lf**2*C_Sf*(g*lr-u[1]*h) + lr**2*C_Sr*(g*lf + u[1]*h))*x[6] \
        #     +mu*m/(I*(lr+lf))*(lr*C_Sr*(g*lf + u[1]*h) - lf*C_Sf*(g*lr - u[1]*h))*x[6] \
        #     +mu*m/(I*(lr+lf))*lf*C_Sf*(g*lr - u[1]*h)*x[2], 
        # (mu/(x[4]**2*(lr+lf))*(C_Sr*(g*lf + u[1]*h)*lr - C_Sf*(g*lr - u[1]*h)*lf)-1)*x[6] \
        #     -mu/(x[4]*(lr+lf))*(C_Sr*(g*lf + u[1]*h) + C_Sf*(g*lr-u[1]*h))*x[6] \
        #     +mu/(x[4]*(lr+lf))*(C_Sf*(g*lr-u[1]*h))*x[2]]
        # 1 / I * (F_yf * math.cos(x[2]) * lf - F_yr * lr),
        # -x[6] + 1 / (m * x[4]) * (F_yf * math.cos(x[2] - x[6]) \
        #     + F_yr * math.cos(x[6]))]
        

    return f


def vehicle_dynamics_st(x, uInit, p, type):
    """
    vehicleDynamics_st - single-track vehicle dynamics
    reference point: center of mass

    Syntax:
        f = vehicleDynamics_st(x,u,p)

    Inputs:
        :param x: vehicle state vector
        :param uInit: vehicle input vector
        :param p: vehicle parameter vector

    Outputs:
        :return f: right-hand side of differential equations
    """

    # set gravity constant
    g = 9.81  # [m/s^2]

    # Retrieve parameters
    mu = p.mu
    lf = p.l_f
    lr = p.l_r
    h = p.h_cg
    m = p.m
    I = p.I_z
    # twf = p.twf
    # twr = p.twr
    # cl_f = p.cl_f
    # cl_r = p.cl_r
    # rho = p.rho
    # A = p.A
    # kroll_f = p.kroll_f
    ax = uInit[1]  # longitudinal acceleration
    delta = uInit[0]  # steering angle
    bank = uInit[2]  # bank angle
    k = uInit[3]

    Bf = p.f_pacejka_b
    Cf = p.f_pacejka_c
    Ef = p.f_pacejka_e
    Br = p.f_pacejka_b
    Cr = p.r_pacejka_c
    Er = p.r_pacejka_e

    # States
    vx = x[3]
    vy = x[4]
    omega = x[5]  # yaw rate
    
    # Compute normal force N with bank angle consideration
    N = m * g * ca.cos(bank) - (m * (vx**2) * k) * ca.sin(bank)

    # Vertical tire forces
    # Fz_f = 0.5 * N * lr / (lf + lr) - 0.5 * h / (lf + lr) * m * ax + 0.25 * cl_f * rho * A * vx ** 2
    # Fz_r = 0.5 * N * lr / (lf + lr) + 0.5 * h / (lf + lr) * m * ax + 0.25 * cl_r * rho * A * vx ** 2
    # gamma_y = vy / vx  # body roll effect

    Fz_f = m * (-ax * h + g * lr) / (lr + lf)
    Fz_r = m * (ax * h + g * lf) / (lr + lf)


    # Fz_fl = Fz_f - kroll_f * gamma_y
    # Fz_fr = Fz_f + kroll_f * gamma_y
    # Fz_rl = Fz_r - (1 - kroll_f) * gamma_y
    # Fz_rr = Fz_r + (1 - kroll_f) * gamma_y

    # Compute lateral tire slip angles
    # a_fl = -ca.atan2((vy + lf * omega) * ca.cos(delta) - vx * ca.sin(delta),
    #                  vx * ca.cos(delta) + (vy + lf * omega) * ca.sin(beta) - 0.5 * twf * omega)
    # a_fr = -ca.atan2((vy + lf * omega) * ca.cos(delta) - vx * ca.sin(delta),
    #                  vx * ca.cos(delta) + (vy + lf * omega) * ca.sin(beta) + 0.5 * twf * omega)
    # a_rl = ca.atan2((lr * omega - vx * ca.sin(beta)),
    #                 (vx * ca.cos(beta) - 0.5 * twr * omega))
    # a_rr = ca.atan2((lr * omega - vx * ca.sin(beta)),
    #                 (vx * ca.cos(beta) + 0.5 * twr * omega))

    # Compute lateral tire forces based on simplified dynamics
    # if type == "pacejka":
    #     F_yf = mu * Fz_fl * ca.sin(p.Cf * ca.arctan(p.Bf * a_fl))
    #     F_fr = mu * Fz_fr * ca.sin(p.Cf * ca.arctan(p.Bf * a_fr))
    #     F_yr = mu * Fz_rl * ca.sin(p.Cr * ca.arctan(p.Br * a_rl))
    #     F_rr = mu * Fz_rr * ca.sin(p.Cr * ca.arctan(p.Br * a_rr))
    # elif type == "linear":
    #     F_yf = mu * Fz_fl * p.C_Sf * a_fl
    #     F_fr = mu * Fz_fr * p.C_Sf * a_fr
    #     F_yr = mu * Fz_rl * p.C_Sr * a_rl
    #     F_rr = mu * Fz_rr * p.C_Sr * a_rr


    a_f = -math.atan((x[4] + x[5] * lf) / x[3]) + delta 
    a_r = -math.atan((x[4] - x[5] * lr) / x[3])

    if type == "pacejka":
        F_yf = mu * Fz_f * ca.sin(Cf * ca.arctan(Bf * a_f))
        F_fr = mu * Fz_f * ca.sin(Cf * ca.arctan(Bf * a_f))
        F_yr = mu * Fz_r * ca.sin(Cr * ca.arctan(Br * a_r))
        F_rr = mu * Fz_r* ca.sin(Cr * ca.arctan(Br * a_r))
    elif type == "linear":
        F_yf = mu * Fz_f * p.C_Sf * a_f
        F_fr = mu * Fz_f * p.C_Sf * a_f
        F_yr = mu * Fz_r * p.C_Sr * a_r
        F_rr = mu * Fz_r * p.C_Sr * a_r

    f = [x[3]*math.cos(x[2]) - x[4]*math.sin(x[2]), 
        x[3]*math.sin(x[2]) + x[4]*math.cos(x[2]), 
        x[5], 
        ax, 
        1/m * (F_yr + F_yf) - x[3] * x[5],
        1/I * (-lr * F_yr + lf * F_yf)] 
        # -mu*m/(x[3]*I*(lr+lf))*(lf**2*C_Sf*(g*lr-u[1]*h) + lr**2*C_Sr*(g*lf + u[1]*h))*x[5] \
        #     +mu*m/(I*(lr+lf))*(lr*C_Sr*(g*lf + u[1]*h) - lf*C_Sf*(g*lr - u[1]*h))*x[6] \
        #     +mu*m/(I*(lr+lf))*lf*C_Sf*(g*lr - u[1]*h)*u[0], 
        # (mu/(x[3]**2*(lr+lf))*(C_Sr*(g*lf + u[1]*h)*lr - C_Sf*(g*lr - u[1]*h)*lf)-1)*x[5] \
        #     -mu/(x[3]*(lr+lf))*(C_Sr*(g*lf + u[1]*h) + C_Sf*(g*lr-u[1]*h))*x[6] \
        #     +mu/(x[3]*(lr+lf))*(C_Sf*(g*lr-u[1]*h))*u[0]]
        # 1 / I * (F_yf * math.cos(u[0]) * lf - F_yr * lr),
        # -x[5] + 1 / (m * x[3]) * (F_yf * math.cos(u[0] - x[6]) \
        #     + F_yr * math.cos(x[6]))]
        

    # # System dynamics (same structure as before)
    # f = [vx * ca.cos(beta + x[4]),
    #      vx * ca.sin(beta + x[4]),
    #      0,  # Placeholder for unused dynamic (e.g., z)
    #      uInit[1],  # Longitudinal acceleration
    #      x[5]]  # Yaw rate

    # # Yaw and lateral dynamics update
    # f.append(1 / I * (F_yf * ca.cos(delta) * lf - F_yr * lr))
    # f.append(-x[5] + 1 / (m * vx) * (F_yf * ca.cos(delta - beta) + F_yr * ca.cos(beta)))

    return f

    #------------- END OF CODE --------------
