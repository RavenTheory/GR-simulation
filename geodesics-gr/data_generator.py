import numpy as np



def generate_initial_satates(num_particles, r_min, r_max):
    """
    Generates random positions and three velocity in spherical coord
    """

    r = np.random.uniform(r_min, r_max, num_particles)
    phi = np.random.uniform(0, 2*np.pi, num_particles)
    theta = np.arccos(np.random.uniform(-1, 1, num_particles))

    spatial_u = np.random.uniform(-0.15, 0.15, (num_particles, 3))
    
    return np.column_stack((r, theta, phi, spatial_u))





states = generate_initial_satates(1000, 3.1, 10.0)


with open("initial_data.csv", "w") as f:
    f.write("r,theta,phi,u1,u2,u3\n")
    for i, state in enumerate(states):
        f.write(",".join(map(str, state)) + "\n")
        #f.write(f"{i}," + ",".join(map(str, state)) + "\n")
