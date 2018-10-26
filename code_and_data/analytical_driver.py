
import numpy as np


def varx_no_norm(t, tau_x):
    return (tau_x / 2) * (1 - np.exp(-t / tau_x) ** 2)


def varx(t, tau_x):
    return (1 - np.exp(-t / tau_x) ** 2)


def cov_rate_YZ(t, gy, gz, tx, ty, tz):
    naut = ((np.exp(-t * (3 / tx + 1 / ty + 1 / tz))) * tx * ty * tz) / (
        ((tx ** 2 - ty ** 2) * (ty + tz) * (tx ** 2 - tz ** 2)) * (-2 * ty * tz + tx * (ty + tz)))

    a = gz ** 2 * (tx + ty) * tz * (
        -2 * np.exp((3 * t) / tx) * (-1 + np.exp(t * (1 / ty + 1 / tz))) * ty ** 2 * tz ** 2 +
        np.exp(t * (1 / tx + 1 / ty + 1 / tz)) * (-1 + np.exp((2 * t) / tx)) * tx ** 3 * (ty + tz) +
        np.exp(t * (1 / tx + 1 / ty)) * (-4 * np.exp(t / tx) + 3 * np.exp(t * (2 / tx + 1 / tz)) + np.exp(
            t / tz)) * tx * ty * tz * (ty + tz) +
        tx ** 2 * (-2 * np.exp((3 * t) / tx) * ty ** 2 + 2 * np.exp(t * (2 / tx + 1 / ty)) * ty ** 2 +
                   np.exp(t * (1 / tx + 1 / ty + 1 / tz)) * ty ** 2 + 4 * np.exp(t * (2 / tx + 1 / ty)) * ty * tz +
                   2 * np.exp(t * (2 / tx + 1 / ty)) * tz ** 2 - np.exp(t * (1 / tx + 1 / ty + 1 / tz)) * tz ** 2 -
                   np.exp(t * (3 / tx + 1 / ty + 1 / tz)) * (ty ** 2 + 4 * ty * tz + tz ** 2)))

    b = gy ** 2 * ty * (tx + tz) * (
        -2 * np.exp((3 * t) / tx) * (-1 + np.exp(t * (1 / ty + 1 / tz))) * ty ** 2 * tz ** 2 +
        np.exp(t * (1 / tx + 1 / ty + 1 / tz)) * (-1 + np.exp((2 * t) / tx)) * tx ** 3 * (ty + tz) +
        np.exp(t * (1 / tx + 1 / tz)) * (-4 * np.exp(t / tx) + 3 * np.exp(t * (2 / tx + 1 / ty)) + np.exp(
            t / ty)) * tx * ty * tz * (ty + tz) +
        tx ** 2 * (2 * np.exp(t * (2 / tx + 1 / tz)) * ty ** 2 -
                   np.exp(t * (1 / tx + 1 / ty + 1 / tz)) * ty ** 2 + 4 * np.exp(t * (2 / tx + 1 / tz)) * ty * tz -
                   2 * np.exp((3 * t) / tx) * tz ** 2 + 2 * np.exp(t * (2 / tx + 1 / tz)) * tz ** 2 +
                   np.exp(t * (1 / tx + 1 / ty + 1 / tz)) * tz ** 2 -
                   np.exp(t * (3 / tx + 1 / ty + 1 / tz)) * (ty ** 2 + 4 * ty * tz + tz ** 2)))

    return naut * (a + b)


def cov_rate(t, g, tx, ty):
    return (g * tx * ty * (
        (1 - np.exp(-2 * t / tx)) * tx + (-1 - np.exp(-2 * t / tx) + 2 * np.exp(-t * (1 / tx + 1 / ty))) * ty)) / (
               tx ** 2 - ty ** 2)


def var_rate(t, g, tx, ty):
    a = (1 / ((tx - ty) ** 2 * (tx + ty))) * np.exp(-t * (1 / tx + 1 / ty))
    b = 4 * g ** 2 * tx ** 2 * ty ** 3
    c = np.exp(t * (-(1 / tx) + 1 / ty)) * g ** 2 * tx ** 2 * ty ** 2 * (tx + ty)
    d = np.exp(t * (1 / tx + 1 / ty)) * (tx - ty) ** 2 * (tx + ty + g ** 2 * tx * ty ** 2)
    e = np.exp(t * (1 / tx - 1 / ty)) * (tx + ty) * (tx ** 2 + ty ** 2 + tx * ty * (-2 + g ** 2 * ty ** 2))

    return a * (b - c + d - e)


def cor2info(corr):
    return -.5 * np.log(1 - corr ** 2)


def varx_no_norm(t, tau_x):
    return (tau_x / 2) * (1 - np.exp(-2 * t / tau_x))


def vary_no_norm(t, tau_x, tau_y, gy):
    u = (np.exp(-2 * t * (1 / tau_x + 1 / tau_y)) * tau_y) / (2 * (tau_x - tau_y) ** 2 * (tau_x + tau_y))

    a = np.exp(2 * t / tau_y) * (-1 + np.exp(2 * t / tau_x)) * gy ** 2 * tau_x ** 4 * tau_y

    b = np.exp(2 * t / tau_x) * (-1 + np.exp(2 * t / tau_y)) * tau_x * tau_y ** 2

    c = np.exp(2 * t / tau_x) * (-1 + np.exp(2 * t / tau_y)) * tau_y ** 3

    d = np.exp(2 * t / tau_x) * (-1 + np.exp(2 * t / tau_y)) * tau_x ** 2 * tau_y * (-1 + gy ** 2 * tau_y ** 2)

    e = tau_x ** 3 * (np.exp(2 * t / tau_x) * (-1 + np.exp(2 * t / tau_y)))

    f = tau_x ** 3 * (np.exp(2 * t / tau_x) - 4 * np.exp(t * (1 / tau_x + 1 / tau_y)) + 2 * np.exp(
        2 * t * (1 / tau_x + 1 / tau_y)) + np.exp(2 * t / tau_y)) * gy ** 2 * tau_y ** 2

    return u * (a - b + c + d + e - f)


def covyz_no_norm(t, tx, ty, tz, gy, gz):
    u = (tx ** 2 * ty * tz) / (
    2 * (tx ** 2 - ty ** 2) * (ty + tz) * (tx ** 2 - tz ** 2) * (-2 * ty * tz + tx * (ty + tz)))

    xnaut = np.exp(-t * (3 / tx + 2 / ty + 1 / tz)) * gy ** 2 * ty * (tx + tz)

    ax = -2 * np.exp(t * (3 / tx + 1 / ty)) * (-1 + np.exp(t * (1 / ty + 1 / tz))) * ty ** 2 * tz ** 2

    bx = np.exp(t * (1 / tx + 2 / ty + 1 / tz)) * (-1 + np.exp(2 * t / tx)) * tx ** 3 * (ty + tz)

    cx = np.exp(t * (1 / tx + 1 / ty + 1 / tz)) * (
    -4 * np.exp(t / tx) + 3 * np.exp(t * (2 / tx + 1 / ty)) + np.exp(t / ty)) * tx * ty * tz * (ty + tz)

    dx_0 = np.exp(t * (1 / tx + 1 / ty)) * tx ** 2
    dx = 2 * np.exp(t * (1 / tx + 1 / tz)) * ty ** 2 - np.exp(t * (1 / ty + 1 / tz)) * ty ** 2 + 4 * np.exp(
        t * (1 / tx + 1 / tz)) * ty * tz - 2 * np.exp(2 * t / tx) * tz ** 2 + 2 * np.exp(
        t * (1 / tx + 1 / tz)) * tz ** 2 + np.exp(t * (1 / ty + 1 / tz)) * tz ** 2 - np.exp(
        t * (2 / tx + 1 / ty + 1 / tz)) * (ty ** 2 + 4 * ty * tz + tz ** 2)

    ynaut = np.exp(-t * (3 / tx + 1 / ty + 2 / tz)) * gz ** 2 * (tx + ty) * tz
    ey = -2 * np.exp(t * (3 / tx + 1 / tz)) * ty ** 2 * (tx ** 2 - tz ** 2)
    fy_0 = np.exp(t * (1 / ty + 1 / tz))
    fy = 2 * np.exp(2 * t / tx) * tx * (ty + tz) * (-2 * ty * tz + tx * (ty + tz))
    fy_1 = 2 * np.exp(t * (2 / tx + 1 / tz)) * (tx - ty) * (
    np.sinh(t / tx) * (-tx * ty * tz + ty * tz ** 2 + tx ** 2 * (ty + tz)) + np.cosh(t / tx) * tz * (
    ty * tz - tx * (2 * ty + tz)))

    return u * (xnaut * (ax + bx + cx + dx_0 * dx) + ynaut * (ey + fy_0 * (fy + fy_1)))


def vary_no_norm(t, tau_x, tau_y, gy):
    u = (np.exp(-2 * t * (1 / tau_x + 1 / tau_y)) * tau_y) / (2 * (tau_x - tau_y) ** 2 * (tau_x + tau_y))

    a = np.exp(2 * t / tau_y) * (-1 + np.exp(2 * t / tau_x)) * gy ** 2 * tau_x ** 4 * tau_y

    b = np.exp(2 * t / tau_x) * (-1 + np.exp(2 * t / tau_y)) * tau_x * tau_y ** 2

    c = np.exp(2 * t / tau_x) * (-1 + np.exp(2 * t / tau_y)) * tau_y ** 3

    d = np.exp(2 * t / tau_x) * (-1 + np.exp(2 * t / tau_y)) * tau_x ** 2 * tau_y * (-1 + gy ** 2 * tau_y ** 2)

    e = tau_x ** 3 * (np.exp(2 * t / tau_x) * (-1 + np.exp(2 * t / tau_y)))

    f = tau_x ** 3 * (np.exp(2 * t / tau_x) - 4 * np.exp(t * (1 / tau_x + 1 / tau_y)) + 2 * np.exp(
        2 * t * (1 / tau_x + 1 / tau_y)) + np.exp(2 * t / tau_y)) * gy ** 2 * tau_y ** 2

    return u * (a - b + c + d + e - f)
if __name__ == "__main__":
    import matplotlib.cm as cm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style('white')
    colors = cm.rainbow(np.linspace(0, 1, 100))
    fig5, ax5 = plt.subplots(2, 2, sharex=True, figsize=(12, 6))
    time_vect = np.linspace(0, 100, 1000)
    gy = 0.1
    gz = 0.1
    tau_x = 5.0
    tau_y = 30.2 / np.log(2)

    tau_z = 30.1 / np.log(2)

    tau_x_range = np.arange(1, 20)[0::2]

    c_m = cm.rainbow

    colors = cm.rainbow(np.linspace(0, 1, 10))
    for z, tau_x in enumerate(tau_x_range):
        sigy = np.sqrt(var_rate(time_vect, gy, tau_x / np.log(2), tau_y))
        sigz = np.sqrt(var_rate(time_vect, gz, tau_x / np.log(2), tau_z))
        covyz = cov_rate_YZ(time_vect, gy, gz, tau_x / np.log(2), tau_y, tau_z)
        corr = covyz / (sigy * sigz)

        ax5[0, 0].plot(time_vect, varx(time_vect, tau_x / np.log(2)), color=colors[z], linewidth=3, alpha=0.75)

        ax5[0, 1].plot(time_vect, cor2info(corr), color=colors[z], linewidth=3, alpha=0.75)

        sigy = np.sqrt(vary_no_norm(time_vect, tau_x / np.log(2), tau_y, gy))
        sigz = np.sqrt(vary_no_norm(time_vect, tau_x / np.log(2), tau_z, gz))
        covyz = covyz_no_norm(time_vect, tau_x / np.log(2), tau_y, tau_z, gy, gz)
        corr = covyz / (sigy * sigz)

        ax5[1, 0].plot(time_vect, varx_no_norm(time_vect, tau_x / np.log(2)), color=colors[z], linewidth=3, alpha=0.75)

        ax5[1, 1].plot(time_vect, cor2info(corr), color=colors[z], linewidth=3, alpha=0.75)

    norm = mpl.colors.Normalize(vmin=min(tau_x_range), vmax=max(tau_x_range))
    c_m = cm.rainbow
    s_m = cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    ax5[1, 0].set_xlabel('Time (Minutes)')
    ax5[0, 1].set_xlabel('Time (Minutes)')
    ax5[0, 0].set_xlabel('Time (Minutes)')
    ax5[1, 1].set_xlabel('Time (Minutes)')

    ax5[0, 0].set_ylabel('Variance')
    ax5[1, 0].set_ylabel('Variance')
    ax5[0, 0].set_title('Diversity of the Activator')
    ax5[1, 0].set_title('Diversity of the Activator')
    ax5[0, 1].set_ylabel('Information (bits)')
    ax5[1, 1].set_ylabel('Information (bits)')

    ax5[0, 1].set_title('Coordination of the Downstream Genes')
    ax5[1, 1].set_title('Coordination of the Downstream Genes')

    # ax5[0].set_xlabel('Time (Minutes)')
    # ax5[0].set_ylabel('Variance')
    # ax5[0].set_title('Activator diversity over time')
    # ax5[1].set_title('Downstream coordination over time')

    fig5.colorbar(s_m, label='Activator Half-Life (minutes)', use_gridspec=False,
                  ax=np.asarray([ax5[0, 0], ax5[0, 1]]).ravel().tolist(), anchor=(2.0, 4.0), orientation='vertical')

    fig5.colorbar(s_m, label='Activator Half-Life (minutes)', use_gridspec=False,
                  ax=np.asarray([ax5[1, 0], ax5[1, 1]]).ravel().tolist(), anchor=(2.0, 4.0), orientation='vertical')

    fig5.tight_layout()
    fig5.subplots_adjust(hspace=.5)
    # fig5.subplots_adjust(wspace=.5)
    # plt.gcf().subplots_adjust(right=0.15)
    fig5.savefig('fig_out/analytical_solutions.pdf', bbox_inches='tight')
