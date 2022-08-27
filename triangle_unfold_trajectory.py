'''

this produces the unfolding plot of a given triangle along a given trajectory

the triangle is defined by its interior angles, and its side of longest length,
its base side, which is placed on the x-axis, its left endpoint at the origin;
the triangle's other two sides face upward, in the +y direction

the trajectory has its foot point somewhere on the base side, and angle
in [0,pi]

in the header section of the main script:
() set the triangle shape: interior angles in fractions of pi (eg [1/10,1/5,
7/10]), and the length of the base side (on x-axis, left edge at origin)
() set the trajectory: trajectory angle (in (0,pi)) and foot point (on the
base side)

'''


import numpy as np
import matplotlib.pyplot as plt
from tri_discontinuity_arcs import Point, TriangleSide, Triangle, plot_sides, \
    two_pt_lin_int, triangle_from_angles


###
# MAIN start
###

# establish triangle shape
largest_side_length = 2.0
angs = np.array([2/10,1/10,7/10])
base_tri = triangle_from_angles(*(np.pi*angs),largest_side_length)

# set trajectory
traj_xvl = 4.331419780219780e-01
traj_ang = 2.233420851148851

# set number of unfoldings
num_unfold = 20

# find the two endpoints of a unit segment along the trajectory and starting at
# foot point on base side
traj_ext = [traj_xvl+np.cos(traj_ang),np.sin(traj_ang)]
traj_pts = [Point(traj_xvl,0.0),Point(*traj_ext)]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

triangle = base_tri
triangle_list = [triangle]
pt_in = None
for jj in range(num_unfold):
    for ix in range(2):
        side = triangle.get_side(ix)
        pt_in = two_pt_lin_int(*side.get_p1p2(),*traj_pts)  # np.array
        if pt_in is not None:
            if side.in_interior(Point(*pt_in)):
                break
    else:
        raise ValueError("problem finding intersection between trajectory and "
                         "triangle sides")
    triangle = triangle.unfold(ix)
    triangle_list.append(triangle)

# plot triangles in infolding
for triangle in triangle_list:
    triangle.plot_triangle()

# plot the trajectory segment
if pt_in is not None:
    plot_sides([TriangleSide(Point(traj_xvl,0.0),Point(*pt_in))])
ax.set_aspect('equal')
plt.axis('off')
plt.show()






