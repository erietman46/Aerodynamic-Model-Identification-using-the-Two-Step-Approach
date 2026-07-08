"""
A simple data file loader
Coen de Visser, TU-Delft, 2026
"""

import scipy.io as sio
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog

# Select datafile to load
Tk().withdraw()
dataname = filedialog.askopenfilename(filetypes=[("MAT files", "*.mat")])

try:
    data = sio.loadmat(dataname)
except Exception:
    raise RuntimeError("Error loading data file")

# Extract variables from the .mat file
# variables in data file:
# t, alpha, Ax, Ay, Az, beta, da, de, dr, dta, dte, dtr, flaps, gamma,
# gear, Mach, p, phi, psi, q, r, Tc1, Tc2, theta, u_n, v_n, vtas, w_n
# string: contains aircraft parameters + variable explanation
var_names = [
    "t", "alpha", "Ax", "Ay", "Az", "beta", "da", "de", "dr",
    "dta", "dte", "dtr", "flaps", "gamma", "gear", "Mach",
    "p", "phi", "psi", "q", "r", "Tc1", "Tc2", "theta",
    "u_n", "v_n", "vtas", "w_n",
]
v = {name: data[name].flatten() for name in var_names}

# --- Plotting ---
plt.close("all")

fig, ax = plt.subplots(num=99)
ax.plot(v["t"], v["p"], "b", label="p [rad/s]")
ax.plot(v["t"], v["q"], "r", label="q [rad/s]")
ax.plot(v["t"], v["r"], "k", label="r [rad/s]")
ax.set_xlabel("time [s]")
ax.legend()

fig, ax = plt.subplots(num=100)
ax.plot(v["t"], v["phi"], "b", label=r"$\phi$ [rad]")
ax.plot(v["t"], v["theta"], "r", label=r"$\theta$ [rad]")
ax.plot(v["t"], v["psi"], "k", label=r"$\psi$ [rad]")
ax.set_xlabel("time [s]")
ax.legend()

fig, ax = plt.subplots(num=101)
ax.plot(v["t"], v["Ax"], "b", label="Ax [m/s²]")
ax.plot(v["t"], v["Ay"], "r", label="Ay [m/s²]")
ax.plot(v["t"], v["Az"], "k", label="Az [m/s²]")
ax.set_xlabel("time [s]")
ax.legend()

fig, ax = plt.subplots(num=102)
ax.plot(v["t"], v["alpha"], "b", label="alpha [rad]")
ax.plot(v["t"], v["beta"], "r", label="beta [rad]")
ax.set_xlabel("time [s]")
ax.legend()

fig, ax = plt.subplots(num=103)
ax.plot(v["t"], v["vtas"], "b", label="VTAS [m/s]")
ax.set_xlabel("time [s]")
ax.legend()

fig, ax = plt.subplots(num=104)
ax.plot(v["t"], v["Tc1"], "b", label="Throttle (left)")
ax.plot(v["t"], v["Tc2"], "r--", label="Throttle (right)")
ax.set_xlabel("time [s]")
ax.legend()

fig, ax = plt.subplots(num=105)
ax.plot(v["t"], v["da"], "b", label="da [rad]")
ax.plot(v["t"], v["de"], "r", label="de [rad]")
ax.plot(v["t"], v["dr"], "k", label="dr [rad]")
ax.set_xlabel("time [s]")
ax.legend()

plt.show()