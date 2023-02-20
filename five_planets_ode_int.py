import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint


# Gravitational Constant

G = 6.67408e-11

# Masses of Sun --> Jupiter

masses = np.array([1988500e24, 0.330e24, 4.87e24,
                   5.97e24, 0.642e24, 1898e24])

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

# Returns dx/dt , dv/dt at different time steps


def derivs(y, t):
    y = y.tolist()
    # Update Position Matrix
    list_ind = 0
    ind2 = 0
    for body in range(N):
        r[body][ind2] = y[list_ind]
        list_ind += 2
    list_ind = 1
    ind2 = 1
    for body in range(N):
        r[body][ind2] = y[list_ind]
        list_ind += 2
    # Update Velocity Matrix
    list_ind = 12
    ind2 = 0
    for body in range(N):
        v[body][ind2] = y[list_ind]
        list_ind += 2
    list_ind = 13
    ind2 = 1
    for body in range(N):
        v[body][ind2] = y[list_ind]
        list_ind += 2

    a = get_acceleration(r, masses)
    for array in a:
        for element in array.tolist():
            y.append(element)
    return y[12:]    # dx/dt and dv/dt values


# Initial Positions and Velocities (Sun --> Jupiter)
r0 = np.array([[0.0, 0.0],
               [0.0, 57.9e9],
               [0.0, 108.2e9],
               [0.0, 149.6e9],
               [0.0, 227.94e9],
               [0.0, 778.6e9]
               ])
v0 = np.array([[0.0, 0.0],
               [47400, 0.0],
               [35000, 0.0],
               [29800, 0.0],
               [24100, 0.0],
               [13100,  0.0]
               ])


y0 = []
for array in r0:
    for element in array.tolist():
        y0.append(element)

for array in v0:
    for element in array.tolist():
        y0.append(element)

# Time Span: Jupiter's orbital period in seconds

t0 = 0
tmax = 377335689.6
dt = 5000
t = np.linspace(t0, tmax, dt)

# Solves the Initial Value Problem
solution = odeint(derivs, y0, t)


x_coords = []
y_coords = []
for coord in range(N*2):
    if coord % 2 == 0:
        x_coords.append(solution[:, coord])
    else:
        y_coords.append(solution[:, coord])


(fig, ax) = plt.subplots(figsize=(4, 4), dpi=190)


planet_colors = ['orange', 'gray', 'darksalmon',
                 'dodgerblue', 'orangered', 'brown']
planet_labels = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter']

# Radii of the sun and the planets

radii = np.array([3.4785e+08, 2.4395e+06, 6.0520e+06,
                  6.3780e+06, 3.3960e+06, 7.1492e+07])

# Marker size in units of Sun Radius

marker_size = [1, 0.00701308,  0.0173983, 0.01833549, 0.00976283, 0.20552537]

# Animation function
# "Text" is the time at a certain frame(seconds) divided by the number of seconds in a day


def animate(frame):
    ax.cla()
    ax.set_xlim(-1e12, 1e12)
    ax.set_ylim(-9e11, 9e11)
    for i in range(N):
        ax.scatter(x_coords[i][frame], y_coords[i][frame], s=75*marker_size[i],
                   c=planet_colors[i], label=planet_labels[i],  marker='o')
        ax.legend(loc='upper left')
        ax.set_facecolor('k')
        ax.text(-0.9e12, -8e11,
                "Time: {:.2f} Days".format(t[frame]/(86400)), color='white')


# Calling the Animating Function(nframe = dt)

anim = animation.FuncAnimation(fig, animate, dt,
                               interval=30, repeat=False)

# Showing a Matplot of the animation
plt.show()

'''
# Saving the Animation; use a writer of your choice
writervideo = animation.FFMpegWriter(fps=60, bitrate=-1)
plt.rcParams['animation.ffmpeg_path'] = r'C:\\Users\\unive\\Downloads\\Programs\\ffmpeg-4.3.1-win64-static\\bin\\ffmpeg.exe'
anim.save('Five Planets.mp4', writer=writervideo)
'''
