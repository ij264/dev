# import matplotlib.pyplot as plt
# import numpy as np
# from scipy.integrate import solve_bvp
#
#
# def rhs(x, y):
#     return np.vstack(
#         (
#             y[1],  # RHS for temperature in upper layer
#             np.repeat(
#                 H_1 / k_1, y[1].size
#             ),  # RHS for temperature derivative in upper layer
#             y[3],  # RHS for temperature in lower layer
#             np.repeat(
#                 H_2 / k_2, y[3].size
#             ),  # RHS for temperature derivative in lower layer
#         )
#     )
#
#
# def bc(ya, yb):
#     return np.array(
#         [
#             yb[0] - T_u,  # BC for temperature at upper layer top
#             ya[2] - T_b,  # BC for temperature at lower layer bottom
#             ya[0] - yb[2],  # BC for continuity of temperature
#             k_1 * ya[1] - k_2 * yb[3],  # BC for continuity of heat flux
#         ]
#     )
#
#
# k_1 = 0.5  # Upper layer thermal conductivity
# k_2 = 1 # Lower layer thermal conductivity
# H_1 = 0  # Upper layer volumetric heating rate
# H_2 = 0  # Lower layer volumetric heating rate
# T_u = 273  # Temperature at upper layer top
# T_b = 1473  # Temperature at lower layer bottom
#
# mesh_node_number = 21
#
# sol = solve_bvp(
#     rhs,
#     bc,
#     np.linspace(0, 1, mesh_node_number),
#     np.vstack(
#         (
#             np.linspace(T_b, 900, mesh_node_number),
#             (T_b - T_u) / 2 * np.ones(mesh_node_number),
#             np.linspace(900, T_u, mesh_node_number),
#             (T_b - T_u) / 2 * np.ones(mesh_node_number),
#         )
#     ),
# )
#
# T_0 = (T_u + T_b) / 2
# print(sol.y[0][0] - T_0)
# plt.plot(
#     np.hstack((sol.y[2], sol.y[0][1:])), np.linspace(20, 0, mesh_node_number * 2 - 1),
#     label=f'k_1 = {k_1}, k_2 = {k_2}'
# )
# plt.axhline(10, linestyle="dotted", linewidth=0.5, color="grey")
# plt.axvline((T_u + T_b)/2, linestyle="dotted", linewidth=0.5, color="grey")
# plt.xlabel("Temperature (K)")
# plt.ylabel("Depth (km)")
# plt.gca().invert_yaxis()
# plt.legend()
# plt.show()

import numpy as np
import matplotlib.pyplot as plt

T_B = 1400
T_T = 300

H = []
k_ratios = []

for i in np.linspace(-2, 2, 101):
    k_ratio = 10 ** i
    d = 7
    z_B = 100

    A_2 = (T_B - T_T) / ((1 / k_ratio - 1) * d + z_B)
    A_1 = A_2 / k_ratio
    B_2 = T_B - A_2 * z_B
    B_1 = T_T


    def T_profile(z):
        if z <= d:
            return A_1 * z + B_1
        else:
            return A_2 * z + B_2


    def heat_flux(z):
        return k_ratio * A_1


    # print(T_profile(7) - 377.0)

    print(f'k_1/k_2 = {k_ratio}. \n The heat flux at the surface is {heat_flux(0):.3} W/m^2')
    H.append(heat_flux(0))
    k_ratios.append(k_ratio)

    z = np.linspace(0, 100, 1000)


plt.plot(k_ratios, H)
plt.xscale('log')
plt.xlabel(r'$k_1/k_2$')
plt.ylabel(r'Surface heat flux $(W/m^2)$')
plt.axvline(1, linestyle="dotted", linewidth=2, color="black")
plt.axhline(11., linestyle="dotted", linewidth=2, color="black")
plt.grid()
plt.show()
# plt.plot([T_profile(z_i) for z_i in z], z, label=f'k_1/k_2 = {k_ratio}')
# plt.ylabel("Depth (km)")
# plt.gca().invert_yaxis()
# plt.xlabel("Temperature (K)")
# plt.axhline(7, linestyle="dotted", linewidth=0.5, color="grey")
# plt.axvline(377, linestyle="dotted", linewidth=0.5, color="grey")
# plt.legend()
# plt.show()
