'''

this routine processes discontinuity boundary vertex arc data from the
tri_discontinuity_arcs module

the output is a list of discontinuity boundary cell polygons and their centroids

given a pickle file with tri_discontinuity_arcs vrt_vis_bas_lst list, ie
a list of 2-tuples, (vertex point,visible interval on base triangle's base
side), this routine converts each 2-tuple to a line segment in the (x,
cotangent(theta)) phase space; the resulting line segments will, within
some degree of precision, always begin and end at another segment (or the
rectangular boundary defining the (x,cot(theta)) phase space window)

the sweep line routine then finds all inter-segment intersections; another
routine connects each of the resulting subsegments that run between
intersection points into polygonal discontinuity boundary cells, where the
triangular billiard bounce sequence is constant at the unfolding level
specified in tri_discontinuity_arcs

three filenames need to be set in the script header:
() disc_arcs_file--this should be a discontinuity arcs file, generated
by tri_discontinuity_arcs
() out_pkl_file--this pickle file will hold the coordinates of segments
composing the polygonal discontinuity boundary cells
() out_cen_file--this text file will hold cell centroid coordinates,
in (x,theta) space

the DEBUG_PLOT flag if set to true will show the polygonal cells as their
cycles are discovered

'''


import pickle
import math
import event_tree  # from associated sweep line repository
from sweep_line import run_sweep_line   # from associated sweep line repository

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib not found; install matplotlib or ammend code for your "
          "plotting modules for full functionality")
    matpltlib_enabled = False
else:
    matpltlib_enabled = True

disc_arcs_file = "" # name of an output file from tri_discontinuity_arcs
out_pkl_file = "" # name of a binary file to write the list of polygon
# segments to
out_cen_file = "" # name of a text file to write the list of polygon
# centroid coordinates to

DEBUG_PLOT = True

TOL_ACC = event_tree.TOL_ACC # tolerance for eg point equality


def angle_compl(ang):
    # helper function; returns pi complement of angle in range [0,2pi)
    return (ang+math.pi)%(2*math.pi)


def centroid(pts_lst):
    # given a list of polygon vertices in order of occurrence around
    # the polygon perimeter, return the centroid;
    # points list is of form [[x1,y1],[x2,y2],...]
    # (https://en.wikipedia.org/wiki/Centroid#Of_a_polygon)
    tp_ls = pts_lst+[pts_lst[0]]
    bs_ls = [tp_ls[ii][0]*tp_ls[ii+1][1]-tp_ls[ii+1][0]*tp_ls[ii][1]
             for ii in range(len(tp_ls)-1)]
    ar = (1/2)*sum(bs_ls)
    if ar==0:   # can happen if eg all points are colinear
        return None
    c_x = (1/(6*ar))*sum(map(lambda x,y:x*y,[tp_ls[ii][0]+tp_ls[ii+1][0]
                                  for ii in range(len(tp_ls)-1)],bs_ls))
    c_y = (1/(6*ar))*sum(map(lambda x,y:x*y,[tp_ls[ii][1]+tp_ls[ii+1][1]
                                  for ii in range(len(tp_ls)-1)],bs_ls))
    return c_x,c_y


def plt_lin(pt1,pt2,clr='green',mrk=False):
    # expects an open plot window
    if not mrk:
        plt.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],color=clr)
    else:
        plt.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],
                 color=clr,marker='o')


def plt_sgs(seg_lst):
    for jj,seg in enumerate(seg_lst):
        plt.plot([seg[0][0],seg[1][0]],[seg[0][1],seg[1][1]],color='blue')
        plt.text(seg[0][0]+(-seg[0][0]+seg[1][0])/2,seg[0][1]+
            (-seg[0][1]+seg[1][1])/2, str(jj),
                 fontsize=9, color='blue')


def plt_cyc(cycle,clr='red'):
    # helper to plot a list of segments, as a polygon
    for sg in cycle:
        plt_lin(sg[0],sg[1],clr=clr)



#@@@
# main
#@@@

with open(disc_arcs_file,'rb') as fp:
    vrt_vis_bas_lst = pickle.load(fp)

# set segment window; this should match the window used in
# tri_discontinuity_arcs to create the arcs pickle file;
# x_lims will be an interval [a,b], a<b, within the base side segment;
# rad_lims will be an angular interval, [c,d], c<d, in radians
x_lims = [] # eg [0.002,1.998]
rad_lims = []   # eg [math.pi/5,2*math.pi/5]

# convert radial limits to limits in the cotangent phase space
y_lims = [1/math.tan(rad_lims[0]),1/math.tan(rad_lims[1])]

# ensure x_lims and y_lims are increasing order
if x_lims[1]<x_lims[0]:
    x_lims = x_lims[::-1]
if y_lims[1]<y_lims[0]:
    y_lims = y_lims[::-1]

# convert [vertex],[x-window] pairs to segments in cotangent space
seg_lst = []    # this will be a list of lists, each sublist a pair of points
# (2-lists) defining segment endpoints; the pair of points always has the
# lowest-x endpoint first, and because of how these segments are ~constructed
# (from vertex,x-window pairs), truly vertical segments are not possible
for ii,tup in enumerate(vrt_vis_bas_lst):
    coa_l = (tup[0][0]-tup[1][0])/tup[0][1]
    coa_r = (tup[0][0]-tup[1][1])/tup[0][1]
    sg_ls = [[tup[1][0],coa_l],[tup[1][1],coa_r]]
    # check for segment fragments--too small to bother with, or wrongly
    # created from arcs routine--
    if (abs(sg_ls[0][0]-sg_ls[1][0])>TOL_ACC and abs(sg_ls[0][1]-sg_ls[1][1])>
            TOL_ACC):
        seg_lst.append(sg_ls)
# add boundary wall segments (L/R, T/B)
sg_ln = len(seg_lst)
lr_wall = [sg_ln,sg_ln+1] # ie [left wall, right wall]
tb_wall = [sg_ln+2,sg_ln+3] # ie [bottom wall, top wall]
for x in x_lims:    # vertical boundary walls
    seg_lst.append([[x,y_lims[0]],[x,y_lims[1]]])
for y in y_lims:    # horizontal boundary walls
    seg_lst.append([[x_lims[0],y],[x_lims[1],y]])


# DEBUG
#print("original seg_lst:")
#print(seg_lst)


# run ~major function for sweep line (this is a fairly substantial function,
# with a lot of code in it)
seg_lst,evs_lst = run_sweep_line(seg_lst)

# obtain angles of segments
ang_dct = {}    # key: value--segment index: angle wrt x axis under left-right
# orientation of the segment, angle in [0,2pi)
for ii,seg in enumerate(seg_lst[:-4]):  # all but boundary segments
    tp_an = math.atan2((seg[1][1]-seg[0][1]),(seg[1][0]-seg[0][0]))
    ang_dct[ii] = tp_an%(2*math.pi)
for sg_ix in lr_wall:   # left/right walls
    ang_dct[sg_ix] = math.pi/2
for sg_ix in tb_wall:   # top/bottom walls
    ang_dct[sg_ix] = 0.0


#@@@
# pre-process segment intersections
#@@@
seg_prc_dct = {ii:event_tree.AVLTree() for ii in range(len(seg_lst))}    #
# dictionary for pre-processing segments; keyed by segment index number
# (from seg_lst)
for ev in evs_lst:  # each event, [x,y],{(sg_ix,L/R/I),()...}
    if len(ev[1])<2:    # a lone segment endpoint should not be allowed here;
        # as well as a point with no segments associated
        print("problem with events in event tree")
        raise ValueError
    for tp in ev[1]:
        # update segment tree:
        dst = math.sqrt((ev[0][0]-seg_lst[tp[0]][0][0])**2 +
                (ev[0][1]-seg_lst[tp[0]][0][1])**2)
        seg_prc_dct[tp[0]].insert(dst,ev[0]) # update each segment's event tree
        # by distance from event point to segment's left endpoint, w/ data the
        # point value [x,y]


#@@@
# process event points for subsegment angular information
#@@@
evt_ang_dct = {}    #
# dictionary for pre-processing; keyed by each event point via (x,y), gives an
# AVL tree with attached segment angles (in [0,2pi)) as values (so, ordered
# CCW); tree entry data is [segment index, next event point ([x,y]) on
# that segment]
for ev in evs_lst:  # each event, [x,y],{(sg_ix,L/R/I),()...}
    evt_ang_dct[(ev[0][0],ev[0][1])] = event_tree.AVLTree()
    for tp in ev[1]:
        sg_ix = tp[0]
        dst = math.sqrt((ev[0][0] - seg_lst[sg_ix][0][0]) ** 2 +
                        (ev[0][1] - seg_lst[sg_ix][0][1]) ** 2)
        if tp[1]=='left' or tp[1]=='internal':
            ang = ang_dct[sg_ix]
            _,ss_pt = seg_prc_dct[sg_ix].hi_neighbor(dst)
            if ss_pt is not None:
                evt_ang_dct[(ev[0][0],ev[0][1])].insert(ang,[sg_ix,ss_pt])
            else:
                print("problem with seg_prc_dct lookup")
                raise ValueError
        if tp[1]=='right' or tp[1]=='internal':
            ang = angle_compl(ang_dct[sg_ix])
            _,ss_pt = seg_prc_dct[sg_ix].lo_neighbor(dst)
            if ss_pt is not None:
                evt_ang_dct[(ev[0][0],ev[0][1])].insert(ang,[sg_ix,ss_pt])
            else:
                print("problem with seg_prc_dct lookup")
                raise ValueError

if DEBUG_PLOT and matpltlib_enabled:
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # plot rectangular wall boundaries
    for jj,seg in enumerate(seg_lst[-4:]):
        plt.plot([seg[0][0],seg[1][0]],[seg[0][1],seg[1][1]],color='blue')
        plt.text(seg[0][0]+(-seg[0][0]+seg[1][0])/2,seg[0][1]+(-seg[0][1]+seg[1][
            1])/2, str(jj+len(seg_lst)-4),
             fontsize=9, color='blue')


#@@@
# prep for polygon construction
#@@@
subsegment_dict = {}    # this will have key:value of type (p1,p2):(p3,p4),
# where everything has been "tuple-ized"--so point pt=[x,y] has become (x,
# y) (to allow dictionary keying)
for sg_ix in range(len(seg_prc_dct)):
    # go through each segment in the tree dictionary, in "both" directions
    sb_sq_ls = [z[1] for z in seg_prc_dct[sg_ix].traverse()]  # defaults to
    # increasing order of [x,y] points along segment w/ index ii; each entry
    # of tree is of form [distance,[x,y]]
    if sb_sq_ls[0] != seg_lst[sg_ix][0]:
        # every segment's left endpoint should occur somewhere in the event
        # list; if not, the first element in sb_sq_ls would not be 0.0
        print("something maybe wrong with subsegments, or event-finding")
        raise ValueError
    if sg_ix != lr_wall[0] and sg_ix != tb_wall[1]:  # don't go CCW off left or
        # top wall in segment L-to-R direction
        for jj in range(len(sb_sq_ls)-1):   # orientation from left endpoint to
            # right endpoint
            angl_tree = evt_ang_dct[tuple(sb_sq_ls[jj+1])]
            ang,dat = angl_tree.lo_neighbor(angle_compl(ang_dct[sg_ix]))
            if ang is None: # next lower angle is across branch cut
                ang,dat = angl_tree.get_max()
            subsegment_dict[(tuple(sb_sq_ls[jj]),tuple(sb_sq_ls[jj+1]))] = \
                (tuple(sb_sq_ls[jj+1]),tuple(dat[1]))
    if sg_ix != lr_wall[1] and sg_ix != tb_wall[0]:   # don't go CCW off
        # right or bottom wall in segment R-to-L direction
        for jj in range(len(sb_sq_ls)-1):   # orientation from left endpoint to
            # right endpoint
            angl_tree = evt_ang_dct[tuple(sb_sq_ls[-jj-2])]
            ang,dat = angl_tree.lo_neighbor(ang_dct[sg_ix])
            if ang is None: # next lower angle is across branch cut
                ang,dat = angl_tree.get_max()
            subsegment_dict[(tuple(sb_sq_ls[-jj-1]),tuple(sb_sq_ls[-jj-2]))] = \
                (tuple(sb_sq_ls[-jj-2]),tuple(dat[1]))

#input("sweep line and subsegment parsing done; press a key: ")


#@@@
# find polygons (closed cycles of linked segments)
#@@@
# cycle_list will have elements of type ((cx,cy),[(p1x,p1y),(p2x,p2y),(p3x,p3y),
# ...]), representing a leading tuple for the centroid, followed by a list of
# the vertex points of the polygon; last and first points make last
# segment
cycle_list = []
while len(subsegment_dict)>0:
    cycle = []
    nxt_seg = next(iter(subsegment_dict))
    cycle.append(nxt_seg)
    while len(cycle) < len(subsegment_dict):
        nxt_seg = subsegment_dict[nxt_seg]
        if nxt_seg == cycle[0]:
            break
        else:
            cycle.append(nxt_seg)
    for ept in cycle:
        del(subsegment_dict[ept])
    ply_pts = [x[0] for x in cycle]
    cycle_list.append((centroid(ply_pts),ply_pts))
    if DEBUG_PLOT and matpltlib_enabled:
        plt_cyc(cycle, clr='purple')
        plt.show()
        plt.pause(0.5)
        #input("press a key: ")

# DEBUG
#print("list of cycles:")
#print(cycle_list)

if DEBUG_PLOT and matpltlib_enabled:
    input("section polygons done; press a key: ")

# save the polygon segment list in a pickle file
with open(out_pkl_file, 'wb') as fp:
    pickle.dump(cycle_list, fp)

# extract list of (x,theta) phase space points for centroids:
phs_cor = [(tup[0][0],math.atan2(1.0,tup[0][1])) for tup in cycle_list] # use
# of atan2 in this way is justified since y values for vertex arc->cotangent
# are always >= 0

# save the polygon centroid coordinates in a text file
with open(out_cen_file, 'w') as fp:
    fp.write("\n".join(["%s, %s" % (x,y) for x,y in phs_cor]))

