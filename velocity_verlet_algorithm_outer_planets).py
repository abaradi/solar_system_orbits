import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint, solve_ivp

# Gravitational Constant

G = 6.67408e-11

# Masses of Sun, Jupiter --> Pluto

masses = np.array([1988500e24, 1898e24, 568e24, 86.8e24, 102e24, 0.0146e24])

# Number of bodies

N = 6  

# Bodies are rows; Coordinates(x and y) are columns

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
    dt = 1000000
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


# Initial positions and velocites of Sun, Jupiter --> Pluto

r0 = np.array([[0.0, 0.0],
               [0.0, 778.6e9],
               [0.0, 1433.5e9],
               [0.0, 2872.5e9],
               [0.0, 4495.1e9],
               [0.0, 5906.4e9]
               ])
v0 = np.array([[0.0, 0.0],
               [13100,  0.0],
               [9700, 0.0],
               [6800, 0.0],
               [5400, 0.0],
               [4700, 0.0]
               ])



# Call the algorithm
solution = integrate(r0, v0, masses)

# Figure Configuration
(fig, ax) = plt.subplots(figsize=(4, 4), dpi=190)

planet_labels = ['Sun', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']

# Radii of the sun and the planets

radii = np.array([3.4785e+08, 7.1492e+07,
                  6.0268e+07, 2.5559e+07, 2.4764e+07, 1.1850e+06])

# Marker size proportional to radii in terms of sun radius

marker_size = [1, 0.20552537,
               0.17325859, 0.07347707, 0.07119161, 0.00340664]

# Animating function 
# "Text" is the time at a certain frame(seconds) divided by the number of seconds in a day

def animate(frame): 
    ax.cla()
    ax.set_xlim(-0.8e13, 0.8e13)
    ax.set_ylim(-0.95e13, 0.95e13)
    for i in range(N):
        ax.scatter(solution[i, 0, frame],
                    solution[i, 1, frame], label=planet_labels[i], s=100*marker_size[i], marker='o')
        ax.legend()
        ax.text(-7e12, -7.5e12,
                "Time: {:.2f} Days".format(frame * 1000000/(86400)), color='k')

# Total number of frames is equal to the total_steps

anim = animation.FuncAnimation(fig, animate, 7825,
                               interval=30, repeat=False)

plt.show()
'''
# Saving the Animation
writervideo = animation.FFMpegWriter(fps=40, bitrate=-1)
plt.rcParams['animation.ffmpeg_path'] = r'C:\\Users\\unive\\Downloads\\Programs\\ffmpeg-4.3.1-win64-static\\bin\\ffmpeg.exe'

anim.save('Outer Planets.mp4', writer=writervideo)
'''
