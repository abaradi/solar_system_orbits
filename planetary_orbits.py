import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint, solve_ivp


# Gravitational Constant

G = 6.67408e-11

# Masses of Sun --> Pluto

masses = np.array([1988500e24, 0.330e24, 4.87e24,
                   5.97e24, 0.642e24, 1898e24, 568e24, 86.8e24, 102e24, 0.0146e24])


# Number of bodies

N = 10

# Bodies in rows and coordinates(x and y) are in columns

r = np.zeros((N, 2))
v = np.zeros((N, 2))
a = np.zeros((N, 2))


# Returns the x and y accleretaions of the planets

def get_acceleration(r, m):
    for i in range(N):
        axj = []
        ayj = []
        for j in range(N):
            if i != j:
                dx = r[i][0] - r[j][0]
                dy = r[i][1] - r[j][1]
                denominator = (dx**2 + dy**2)**(1.5)
                axj.append(-G * ((m[j] * dx) / denominator))
                ayj.append(-G * ((m[j] * dy) / denominator))
                a[i][0] = sum(axj)
                a[i][1] = sum(ayj)
    return a

# Velocity Verlet Algorithm updates the position matrices

def integrate(r, v, m):
    t = 0
    tmax = 7.82438e+9
    dt = 100000
    total_steps = int(tmax/dt)
    # Initial Acceleration
    a = get_acceleration(r, m)
    # Planets on axis 0, Coordinates(x and y) on axis 1, Number of Time Step on axis 3
    solution = np.zeros((N, 2, total_steps + 1))
    solution[:, :, 0] = r
    for i in range(total_steps):
        v = v + a * (dt / 2)
        r = r + v * (dt)
        a = get_acceleration(r, m)
        v = v + a * (dt / 2)
        t = t + dt
        solution[:, :, i + 1] = r
    return solution


# Initial Positions and Velocities (Sun --> Pluto)
r0 = np.array([[0.0, 0.0],
               [0.0, 57.9e9],
               [0.0, 108.2e9],
               [0.0, 149.6e9],
               [0.0, 227.9e9],
               [0.0, 778.6e9],
               [0.0, 1433.5e9],
               [0.0, 2872.5e9],
               [0.0, 4495.1e9],
               [0.0, 5906.4e9]
               ])
v0 = np.array([[0.0, 0.0],
               [47.4e3, 0.0],
               [35.0e3, 0.0],
               [29.8e3, 0.0],
               [24.1e3, 0.0],
               [13100,  0.0],
               [9700, 0.0],
               [6800, 0.0],
               [5400, 0.0],
               [4700, 0.0]
               ])


# Call the algorithm
solution = integrate(r0, v0, masses)


(fig, ax) = plt.subplots(figsize=(8.5, 8.5))

# Figure Customization

planet_colors = ['orange', 'green', 'blue', 'red',
                 'yellow', 'brown', 'gold', 'turquoise', 'lightskyblue', 'blue']
planet_labels = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars',
                 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']

for i in range(N):
    if i == 0:
        ax.scatter(solution[i, 0, 0], solution[i, 0, 0])
    ax.plot(solution[i, 0, 1:], solution[i, 1, 1:],
                    c=planet_colors[i], label=planet_labels[i], linewidth=0.5)
    #ax.legend(loc='upper left')
#        ax.set_facecolor('k')

# Showing a Matplot of the animation
#plt.show()
fig.savefig("Planetary Orbits.png")
