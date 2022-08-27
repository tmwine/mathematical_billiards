'''

generates discontinuity boundary arcs for triangle unfolding at fixed depth

this allows inputting of a triangle's interior angles, a base (longest)
side length, a phase space window in (x,theta) phase space, and the
number of triangle unfoldings, and returns (and potentially plots) the
associated discontinuity boundary cells (see for instance A. Katok's,
The Growth Rate for the Number of Singular and Periodic Orbits for a
Polygonal Billiard, Comm. Math. Phys. 111(1): 151-160 (1987))

to set up, in the header area of "main" script portion:
() input the angle list (eg for a triangle with angles pi/10,pi/5,7*pi/10,
enter [1/10,2/10,7/10])
() set the base side length (this is arbitrary and has no impact on output; a
typical value is 2.0)
() set the zoom window:
    base_subinterval = [left_endpoint,right_endpoint], where the
endpoints must be within the base side's limits
    radial_target = [CW lower limit, CCW upper limit]
    (note, the zoom window defaults to the full base_subinterval,
and radial_target [0,pi])
() set the iteration depth, num_unfold--this is how many unfolding levels to do

if wanting to export the unfolded vertex arcs data, then set to an appropriate
filename for out_pkl_fle, the pickle file for storing vrt_vis_bas_lst

there is also an optional plotting section at the end of "main"; if only
wanting to export the vertex arcs data, this, and the matplotlib dependency
is not necessary

'''


import numpy as np
import matplotlib.pyplot as plt
import copy
import pickle
import math
import time

INT_EPS = 1.0e-6 # reference accuracy for approximating the interior
# of an interval; note this cannot be very small, otherwise small angle
# effects get lost to precision
TOL_ACC = 1e-12 # general tolerance used in some functions
TOL_TOO = max(1e-2*INT_EPS,1e-12)  # another tolerance, used eg in checking
# point on interior of side segment; note this should be a few orders of
# magnitude smaller than INT_EPS


class Point:

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.xy_lst = [x,y]

    def __eq__(self,other):
        if math.isclose(self.get_x(),other.get_x(),rel_tol=0.0,
         abs_tol=TOL_ACC) and math.isclose(self.get_y(),other.get_y(),
         rel_tol=0.0,abs_tol=TOL_ACC):
            return True
        else:
            return False

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_xy(self):
        # returns a list, [x,y]
        return self.xy_lst

    def get_np_arr(self):
        return np.array([self.x,self.y])


class RadialInterval:
    '''
    this accepts two angles that define a radial interval;
    the radial interval will go from the first angle, then in a CCW direction
    to the second;
    angles will be internally stored in branch [0,2pi]
    eg pi/4,pi/2 will include everything between 45 and 90 degrees
    eg 3pi/2,pi/4 will include everything from the negative y axis,
    CCW around to the +45 degree line
    eg pi/2,pi/4 would assume the "long way around," for a total angular
    sweep of 7pi/4
    '''

    def __init__(self,theta_1,theta_2):
        # lower and upper are defined wrt CCW direction; so theta_1 starts,
        # go CCW, and theta_2 ends; 0 is the angle along positive x axis
        self.lower_angle = np.mod(theta_1,2*np.pi)
        self.upper_angle = np.mod(theta_2,2*np.pi)
        self.dub_la = None
        self.dub_ua = None
        self.double()

    def double(self):
        # 'doubles' the radial interval, by adding pi to the whole interval;
        # this treats the interval like lines, vs rays
        self.dub_la = np.mod(np.pi+self.lower_angle,2*np.pi)
        self.dub_ua = np.mod(np.pi+self.upper_angle,2*np.pi)

    def get_angs(self):
        return [self.lower_angle,self.upper_angle]  # 1st element in list
        # always the 'lower' (more clockwise) angle

    def get_dubs(self):
        return [self.dub_la,self.dub_ua] # 1st element in list
        # always the 'lower' (more clockwise) angle

    def check_inside(self,arc_pos):
        # checks if value arc_pos is inside the radial interval; auto-mods
        # arc_pos into [0,2*pi] interval
        # note, this includes the closed interval endpoints--so arc_pos can
        # be on the radial interval boundary and be considered "in" the interval
        arc_pos = np.mod(arc_pos,2*np.pi)
        if self.lower_angle > self.upper_angle:   # angle straddles 0
            return arc_pos>=self.lower_angle or arc_pos<=self.upper_angle
        else:
            return self.lower_angle<=arc_pos and arc_pos<=self.upper_angle

    def get_line_intervals(self):
        # for a given radial interval, extend the rays determining the radial
        # interval as lines; consider the union of these (now two/paired)
        # radial intervals (roated 180 deg from each other); these amount to
        # a line sweep interval, which can be completely defined by one or
        # two line sweeps (intervals) in [0,pi);
        # intervals are returned in increasing order;
        # this returns lists (of lists) (not numpy arrays)--sublists contain
        # floats in [0,pi]
        tmp_upper = self.upper_angle
        if self.upper_angle<self.lower_angle:
            tmp_upper = self.upper_angle + 2*np.pi
        if tmp_upper-self.lower_angle >= np.pi:
            # if the radial interval is >= pi, line sweep is all of
            # S^1; so return the whole [0,pi] interval
            return [[0.0,np.pi]]
        if self.lower_angle >= np.pi: # must be in [pi,2*pi)
            if not self.upper_angle>=self.lower_angle: # upper must be in [0,
                # lower_angle-pi); straddles 0
                if self.upper_angle==0.0:
                    return [[self.lower_angle-np.pi,np.pi]]
                else:
                    return [[0.0,self.upper_angle],[self.lower_angle-np.pi,
                                                   np.pi]]
            else:   # lower and upper are in [pi,2*pi)
                return [[self.dub_la,self.dub_ua]]
        else:   # lower is in [0,pi)
            if self.upper_angle>=np.pi: # straddles pi
                if self.upper_angle==np.pi:
                    return [[self.lower_angle,np.pi]]
                else:
                    return [[0.0,self.upper_angle-np.pi],[self.lower_angle,
                                                          np.pi]]
            else:
                return [[self.lower_angle,self.upper_angle]]


class TriangleSide:

    def __init__(self, p1, p2):   # p's are point objects
        self.p1 = p1
        self.p2 = p2
        self.p1_int = Point(0.0,0.0)
        self.p2_int = Point(0.0,0.0)
        self.int_tol = 1e-5 # internal accuracy level for saving time re
        # checking point in segment interior
        self.set_interior_lims()

    def swap_pts(self):
        # swaps order of pts
        pass

    def get_p1(self):
        return self.p1

    def get_p2(self):
        return self.p2

    def get_p1p2(self):
        return [self.p1,self.p2]

    def get_p1_int(self):
        return self.p1_int

    def get_p2_int(self):
        return self.p2_int

    def get_p1p2_int(self):
        return [self.p1_int,self.p2_int]

    def get_midpt(self):
        # return Point object at midpoint of the segment
        int_vec = self.p2.get_np_arr()-self.p1.get_np_arr()
        return Point(*(self.p1.get_np_arr()+int_vec/2))

    def in_interior(self,pt):
        # determines if Point object pt is "on" the segment; determines to
        # some acceptable tolerance; effectively checks for proximity to a
        # ~tolerance band around the segment
        if self.p1==pt or self.p2==pt:
            return True
        v1 = make_vec(self.p1,pt,norm=False)
        nv1 = np.linalg.norm(v1)
        v1u = v1/nv1
        vin = make_vec(self.p1,self.p2,norm=False)
        nvi = np.linalg.norm(vin)
        viu = vin/nvi
        dt = np.dot(v1u,viu)
        #if abs(dt) < 1.0-self.int_tol:  # if angle sufficiently far from 0 or
            # pi, return
        if abs(dt) < 1.0-self.int_tol:
            return False
        sg_ps = nv1*dt  # (scalar) position along segment of the projection
        if sg_ps-nvi > TOL_TOO or sg_ps < -TOL_TOO:
            return False
        pd = np.linalg.norm(self.p1.get_np_arr()+sg_ps*viu - pt.get_np_arr())
        if pd > TOL_TOO:
            return False
        return True

    def set_interior_lims(self):
        # determines points just inside endpoints; simulate open interval
        vec_p1p2 = make_vec(self.p1,self.p2,norm=True) # from p1 to p2
        self.p1_int = Point(self.p1.get_x()+INT_EPS*vec_p1p2[0],
                       self.p1.get_y()+INT_EPS*vec_p1p2[1])
        self.p2_int = Point(self.p2.get_x()-INT_EPS*vec_p1p2[0],
                       self.p2.get_y()-INT_EPS*vec_p1p2[1])


class Triangle:

    '''
    expects TriangleSide objects
    '''

    def __init__(self, side1, s_num1, side2, s_num2, base_side=2):
        # holds sides in list, in order [side1,side2,side3];
        # holds vertices in list, [side3<>side1 vertex, side1<>side2 vertex,
        # side2<>side3 vertex];
        # note sides are indexed internally from 0;
        # s_num- are which symbol to assign to resp sides; s_num is an index
        # in {0,1,2}, corresponding to 'a','b','c'
        tmp_st = {0,1,2}-{s_num1,s_num2}
        if len(tmp_st) != 1:
            print("problem with Triangle s_num indices")
            raise ValueError
        # sym_nums will assign symbol number sym_nums[0] to first side in
        # self.sides, [1] to second side, etc.:
        self.sym_nums = [s_num1,s_num2,list(tmp_st)[0]]
        s1_ps = side1.get_p1p2()
        s2_ps = side2.get_p1p2()
        cv_ct = 0
        cv_ix = (-1,-1)
        for ii in range(2):
            for jj in range(2):
                if s1_ps[ii]==s2_ps[jj]:
                    cv_ix = (ii,jj) # note prospective common vertex,
                    # in terms of side point indices
                    cv_ct += 1
        if cv_ct != 1:
            print("bad constructor call to Triangle object")
            raise ValueError
        self.vertices = []  # triangle vertices as point objects
        self.vertices.append(s1_ps[1-cv_ix[0]]) # first vertex in list is
        # the point in side1 not at the side1<->side2 vertex
        self.vertices.append(s1_ps[cv_ix[0]])   # next vertex in list is the
        # found common vertex between side1 and side2 (ie the other point on
        # side1)
        self.vertices.append(s2_ps[1-cv_ix[1]]) # last vertex is the
        # remaining point in side2
        self.sides = []
        self.sides.append(side1)
        self.sides.append(side2)
        self.sides.append(TriangleSide(self.vertices[2],self.vertices[0]))
        self.base_side = base_side  # this may not be used; base is typically
        # always side indexed by 2 (starting w/ index 0)

    def __eq__(self, other):
        ss = {0,1,2}
        so = {0,1,2}
        while ss and so:
            el1 = next(iter(ss))
            for el2 in so:
                if self.get_vertex(el1)==other.get_vertex(el2):
                    break
            else:
                return False
            ss.remove(el1)
            so.remove(el2)
        return True

    def change_base(self,ix):
        if ix in range(3):
            self.base_side = ix
        else:
            print("bad value passed to change_base in Triangle object")
            raise ValueError

    def unfold(self,ix):
        # returns a new triangle object, based on unfolding along side "ix";
        # this always makes the two new sides the 1st and 2nd side of the
        # newly created triangle object (ie sides internally indexed 0 and 1)
        if ix not in range(3):
            print("bad value passed to unfold in Triangle object")
            raise ValueError
        sym_vec = make_vec(self.vertices[ix],self.vertices[(ix+1)%3],norm=True)
        # sym_vec is a unit vector along line of symmetry to do reflection /
        # unfolding over; tail is at vertex ix
        osi = (ix+1)%3  # side index for other side, one index higher mod 3;
        # this side will be of points self.vertices[osi] to self.vertices[(
        # osi+1) mod 3]; tail is at vertex on line of symmetry
        sid_vec = make_vec(self.vertices[osi],self.vertices[(osi+1)%3],
                           norm=False)
        sym_lin_ref = np.dot(sid_vec,sym_vec)*sym_vec # vector representing
        # projection of sid_vec onto line of symmetry; referenced to
        # vertices[osi]
        sv_ccw = np.array([-sym_vec[1], sym_vec[0]])  # 90 degree CCW
        # rotation of unit vector sym_vec / perpendicular
        prp_lin_ref = np.dot(sid_vec,sv_ccw)*sv_ccw  # vector representing
        # projection of sid_vec onto perpendicular to line of symmetry;
        # referenced to vertices[osi]
        new_vrt = self.vertices[osi].get_np_arr()+sym_lin_ref-prp_lin_ref
        new_pt = Point(*new_vrt)
        new_side1 = TriangleSide(self.vertices[ix],new_pt)
        new_side2 = TriangleSide(self.vertices[osi],new_pt)
        return Triangle(new_side1,self.sym_nums[(ix-1)%3],new_side2,
                        self.sym_nums[osi])

    def get_side(self,ix):
        if ix not in range(3):
            print("bad value passed to Triangle get_side")
            raise ValueError
        return self.sides[ix]

    def get_side_nums(self):
        return self.sym_nums

    def get_vertex(self,ix):
        # returns point object
        # note, ix=1 (out of {0,1,2}) is the vertex between side1 and side2
        # (ie sides indexed 0 and 1)
        return self.vertices[ix]

    def vertex_dump(self):
        # returns all 3 vertices as a 3x2 numpy array; vertices (rows) are
        # ordered by side-labelling in self.sym_nums;
        # recall self.sym_nums[0] corresponds to side self.sides[0], etc
        # so eg sym_nums=[2,0,1] means sides[0] is assigned syn num 2, etc.
        out_arr = np.zeros([3,2])
        sds_ixs = [self.sym_nums.index(ii) for ii in range(3)] # inverts
        # sym_nums, so eg sym_nums=[2,0,1] becomes [1,2,0] (so in terms of
        # sides indices to symbol #--side index 1 has symbol # 0 (0th spot in
        # list), etc); this allows converting from symbol numbers to actual
        # triangle vertices
        tst_arr = [[[2,0],[0,2]], [[0,1],[1,0]], [[1,2],[2,1]]]
        for ii in range(3):
            adj_sds = [sds_ixs[ii],sds_ixs[(ii+1)%3]]   # "true" indices of
            # sides from symbol # ii to symbol # ii+1 (mod 3):
            # [2,0] or [0,2] gets vertices[0]
            # [0,1] or [1,0] gets vertices[1]
            # [1,2] or [2,1] gets vertices[2]
            for jj,ls in enumerate(tst_arr):
                if adj_sds==ls[0] or adj_sds==ls[1]:
                    break
            else:
                print("error in triangle vertices dump")
                raise ValueError
            out_arr[(ii+1)%3,:] = self.get_vertex(jj).get_np_arr()
        return out_arr

    def plot_triangle(self,clr='blue'):
        for ii in range(3):
            plt_lin(self.vertices[ii],self.vertices[(ii+1)%3],clr=clr)


def plt_lin(pt1,pt2,clr='blue',mrk=False):
    # expects an open plot window
    if not mrk:
        plt.plot([pt1.get_x(),pt2.get_x()],[pt1.get_y(),pt2.get_y()],color=clr)
    else:
        plt.plot([pt1.get_x(), pt2.get_x()], [pt1.get_y(), pt2.get_y()],
                 color=clr,marker='o')


def pp_dist(pt1,pt2):
    # determine 2-norm distance between two point objects
    return np.sqrt((pt1.get_x()-pt2.get_x())**2+(pt1.get_y()-pt2.get_y())**2)


def make_vec(pt1,pt2,norm=False):
    # returns the vector (2-element np array) with tail at pt1 and head at pt2
    if pt1==pt2:
        return None
    ret_vec = pt2.get_np_arr()-pt1.get_np_arr()
    if norm:
        return ret_vec/np.linalg.norm(ret_vec)
    else:
        return ret_vec


def triangle_from_angles(ang1,ang2,ang3,lsl):

    # constructs a triangle object, as 2 "outward facing" side objects,
    # w/ longest side, of length lsl, on the x axis, w/ leftmost
    # vertex at origin;
    # side lengths decrease in counterclockwise order

    ang_arr = np.array([ang1,ang2,ang3])
    if abs(sum(ang_arr)-np.pi)>0.001:
        print("problem with angles")
        return None
    ang_arr = np.sort(ang_arr)  # increasing order
    h2 = lsl / (
                np.sin(ang_arr[0]) * np.cos(ang_arr[1]) / np.sin(ang_arr[1])
         + np.cos(ang_arr[0]))
    h1 = (np.sin(ang_arr[0]) / np.sin(ang_arr[1])) * h2

    vrt = np.array([h1 * np.cos(ang_arr[1]), h1 * np.sin(ang_arr[1])])
    llv = Point(0.0,0.0) # lower left vertex (at origin)
    lrv = Point(lsl,0.0) # lower right vertex
    upv = Point(*vrt)   # upper vertex (what was just computed)
    ul_side = TriangleSide(llv,upv) # should be shortest side; give symbol 2
    ur_side = TriangleSide(upv,lrv) # middle-length side; give symbol 1
    return Triangle(ul_side,2,ur_side,1)


def pt_line_rad(pt,seg):
    # receives a point and a segment, and returns a radial interval,
    # being the S^1 arc occupied by the segment from the perspective at the
    # point;
    # it's possible pt is a segment end point, or (pt) lies within the
    # segment--will return None in this case;
    # note, this does not do radial angle "doubling" (that's internal to radial
    # arc object)

    if seg.in_interior(pt):
        return None # can't really get solid angle, since pt is effectively
        # on the line segment
    # normalized vectors from pt to segment end points:
    vp1 = make_vec(pt,seg.get_p1(),norm=True)
    vp2 = make_vec(pt,seg.get_p2(),norm=True)
    tmp = np.dot(vp1,vp2)
    if abs(tmp) > 1.0:  # vp1 lies along x axis
        tmp = round(tmp)
    ang_dff = np.arccos(tmp)    # always in [0,pi]
    tmp = vp1[0]
    if abs(tmp) > 1.0:  # vp1 lies along x axis
        tmp = round(tmp)
    an1 = np.arccos(tmp)  # angle of vp1
    if vp1[1]<0.0:  # lies in -y halfplane
        an1 = -an1
    v1ccw = np.array([-vp1[1],vp1[0]])  # 90 degree CCW rotation of vp1
    # orientation vector for segment
    vsg = make_vec(seg.get_p1(),seg.get_p2(),norm=True)
    tmp = np.dot(v1ccw,vsg) # positive if segment end p2 lies CCW of vp1 radial

    # DEBUG
    #if np.isnan(an1) or np.isnan(ang_dff):
    #    pass

    if tmp > 0.0:
        return RadialInterval(an1,an1+ang_dff)
    else:
        return RadialInterval(an1-ang_dff,an1)


def fetch_edges(st_ix,edges_lst):
    # fetch edge history for the element in edges_lst[-1][st_ix]; return a
    # list of edges; note, order here may be important--this will return
    # edges in child-parent order--so the first edge in out_lst should be the
    # latest unfolding, while the last edge in out_lst (out_lst[-1]) should
    # be the "base" edge of the starting triangle / polygon
    out_lst = []
    ct = 1
    ix = st_ix
    #print("fetch edges:",end=" ")   # DEBUG
    while ct <= len(edges_lst):
        #print(edges_lst[-ct][ix][1],end="; ")  # DEBUG; symbol number
        out_lst.append(edges_lst[-ct][ix][2])
        ix = edges_lst[-ct][ix][0]
        ct += 1
    #print() # DEBUG
    return out_lst


def arc_intersect(arc_1,arc_2):
    # finds intersection of radial intervals arc_1 and arc_2; returns a list
    # of radial intervals, None if DNE; if one or more radial intervals are
    # > 180 degrees, this may return multiple intersections;
    # angs_. will be [lower_angle,upper_angle] w/ lower_angle the more
    # clockwise angle; angles are in [0,2*pi]
    angs_1 = arc_1.get_angs()
    angs_2 = arc_2.get_angs()
    # check if arcs are identical
    if np.allclose(angs_1,angs_2,rtol=0.0,atol=TOL_ACC):
        return [RadialInterval(*angs_1)]
    chk_2_lo_in_1 = arc_1.check_inside(angs_2[0])   # "inside" includes on
    # the (closed) boundary of the radial interval
    chk_2_hi_in_1 = arc_1.check_inside(angs_2[1])
    chk_1_lo_in_2 = arc_2.check_inside(angs_1[0])
    chk_1_hi_in_2 = arc_2.check_inside(angs_1[1])
    if not (chk_2_lo_in_1 or chk_2_hi_in_1 or chk_1_lo_in_2 or chk_1_hi_in_2):
        # neither interval is within the other; no overlap
        return None
    chk_1_in_2 = chk_1_lo_in_2 and chk_1_hi_in_2
    chk_2_in_1 = chk_2_lo_in_1 and chk_2_hi_in_1
    if chk_1_in_2 and chk_2_in_1:
        # both intervals are within each other; can happen if at least
        # one interval is > 180 degrees
        ret_lst = []
        if abs(angs_1[0]-angs_2[1])>TOL_ACC:
            ret_lst.append(RadialInterval(angs_1[0],angs_2[1]))
        if abs(angs_2[0]-angs_1[1])>TOL_ACC:
            ret_lst.append(RadialInterval(angs_2[0],angs_1[1]))
        if ret_lst:
            return ret_lst
        else:
            print("problem in arc_intersect")
            raise ValueError
    elif chk_1_in_2:
        return [RadialInterval(angs_1[0],angs_1[1])]
    elif chk_2_in_1:
        return [RadialInterval(angs_2[0],angs_2[1])]
    # if here, there is some overlap, but neither interval is completely
    # contained in the other
    if chk_1_lo_in_2: # (which by here must mean chk_1_hi_in_2 is False,
        # and that chk_2_hi_in_1 is True)
        return [RadialInterval(angs_1[0],angs_2[1])]
    elif chk_2_lo_in_1: # (by here it must be that chk_2_hi_in_1 is False,
        # and that chk_1_hi_in_2 is True)
        return [RadialInterval(angs_2[0],angs_1[1])]
    print("problem in arc_intersect")
    raise ValueError    # if made it to here, something's wrong


def get_vertex_vis(vt_pt,sde_lst):
    # receives a point object (corresponding to triangle vertex), and a list
    # of side objects in order of unfolding (0th index is latest unfolding,
    # and sde_lst[-1] is the originating triangle's base);
    # returns a list of arcs (normally, a single element) comprising all solid
    # angles occupied by all segments in sde_lst visible from vt_pt (ie solid
    # angle intersection)

    if sde_lst:
        tmp_arc = pt_line_rad(vt_pt,sde_lst[0]) # initialize
        if tmp_arc is None:
            return None
        else:
            rad_arcs = [tmp_arc] # initialize
    else:
        return None
    for seg in sde_lst[1:]:
        plr = pt_line_rad(vt_pt,seg)
        if plr is None:
            return None # vt_pt is on some segment seg
        arc_lst = []    # intersections of all arcs in rad_arcs with plr arc
        for rad_arc in rad_arcs:
            tmp_arc = arc_intersect(rad_arc,plr)
            if tmp_arc:
                arc_lst += tmp_arc
        if not arc_lst: # at any point if intersection is empty, return
            return None
        rad_arcs = arc_lst
    return rad_arcs


def vtx_arc_int_con(vt_pt,arc,base_subint,ang_win):
    # receives a point object, corresponding to triangle vertex point,
    # and a radial interval object, corresponding to angular visible window
    # from that vertex point through all unfolded sides of relevance to this
    # point, and a base subinterval TriangleSide object (usually along x-axis),
    # and an angular limiting window as a radial interval object (the angles
    # the endpoints of the vertex make with respect to the interval endpoints
    # must be within ang_win);
    # returns a 2-list, an interval within base_subint that meets both
    # angular conditions with respect to the vertex;
    # it's possible the angular window prevents any sightlines--then returns
    # None
    base_side = copy.deepcopy(base_subint)
    tmp_p1,tmp_p2 = base_side.get_p1p2()
    base_left = tmp_p1.get_x()
    base_right = tmp_p2.get_x()
    if base_right < base_left:
        tmp_bse = base_left
        base_left = base_right
        base_right = tmp_bse
    base_arc = pt_line_rad(vt_pt,base_side) # this really shouldn't need
    # checking (right? arc was already checked re sightline to whole base
    # segment base_subint)
    if base_arc is None: # unlikely, but if vt_pt is within base_side
        return None
    com_arc = arc_intersect(arc,base_arc)[0]
    # com_arc will be a radial interval "within" the vertex angle of the
    # triangle formed between point vt_pt and the base segment on x-axis
    if ang_win is not None:
        tmp_ins = arc_intersect(RadialInterval(*ang_win.get_dubs()),com_arc)
        if tmp_ins is None:
            return None
        else:
            win_arc = tmp_ins[0] #
    # refine com_arc to fit within window arc limits
    else:
        win_arc = com_arc
    arc_angs = win_arc.get_angs()
    ang_pt1 = Point(vt_pt.get_x()+np.cos(arc_angs[0]),vt_pt.get_y()+ np.sin(
        arc_angs[0]))
    ang_pt2 = Point(vt_pt.get_x()+np.cos(arc_angs[1]),vt_pt.get_y()+ np.sin(
        arc_angs[1]))
    in_st = two_pt_lin_int(vt_pt,ang_pt1,*base_side.get_p1p2()) # start of
    # interval, as numpy 2-element array (x,y)
    in_nd = two_pt_lin_int(vt_pt,ang_pt2,*base_side.get_p1p2()) # end of
    # interval, as numpy 2-element array (x,y)
    x_0 = in_st[0]
    x_1 = in_nd[0]
    if abs(in_st[1])>TOL_ACC or abs(in_nd[1])>TOL_ACC or x_0>x_1: #
        # should have y~=0.0, and x of in_st should be < x of in_nd
        raise ValueError
    # for any precision issues:
    if x_0<base_left:
        x_0=base_left
    if x_1>base_right:
        x_1 = base_right
    return [x_0,x_1]


def two_pt_lin_int(pt_1_1,pt_1_2,pt_2_1,pt_2_2):
    # finds the intersection (if any) between the two lines, each line
    # defined by two points on it;
    # returns numpy array representing the point of intersection
    pt_lst = [[pt_1_1,pt_1_2],[pt_2_1,pt_2_2]]
    aa = np.zeros([2,2])
    bb = np.zeros(2)
    for ii in range(2): # loop through each line
        two_pts = pt_lst[ii]
        aa[ii,0] = two_pts[0].get_y()-two_pts[1].get_y()    # x coeff
        aa[ii,1] = two_pts[1].get_x()-two_pts[0].get_x()    # y coeff
        bb[ii] = aa[ii,1]*two_pts[1].get_y()+aa[ii,0]*two_pts[1].get_x()
    try:
        ret_arr = np.linalg.solve(aa,bb)
    except np.linalg.LinAlgError as err:
        if 'Singular matrix' in str(err):
            ret_arr = None
        else:
            raise
    return ret_arr


def side_check_sor(sde,ray_ang,x_cor,l_or_r):
    # check if side "sde" is to the l_or_r side of a ray w/ angle ray_ang
    # emanating from point (x_cor,0.0); sde also must have both points w/ y
    # coordinate >= 0;
    # l_or_r is either 'right' or 'left';
    # ray_ang is expected in [0,pi]
    pts = sde.get_p1p2()
    y_vals = [pts[0].get_y(),pts[1].get_y()]
    if y_vals[0] < 0.0 or y_vals[1] < 0.0:
        return False
    if ray_ang == np.pi:
        if l_or_r == 'right':
            return True # everything is to the "right" of -x dir ray
        else:
            return False # nothing is to the "left" of -x dir ray
    elif ray_ang == 0.0:
        if l_or_r == 'right':
            return False # nothing is to the "right" of +x dir ray
        else:
            return True # everything is to the "left" of +x dir ray
    x_vals = [pts[0].get_x(),pts[1].get_x()]
    for ii in range(2):
        x_dff = y_vals[ii]/np.tan(ray_ang)
        if l_or_r == 'right':
            if x_vals[ii] < x_dff+x_cor:
                break
        else:
            if x_dff+x_cor < x_vals[ii]:
                break
    else:
        return True
    return False


def hull_init(in_seg):
    # initialize hulls w/ base segment endpoints;
    # this expects in_seg to be a base segment on the x axis
    l_hull = []
    r_hull = []
    ptl = in_seg.get_p1().get_xy()
    ptr = in_seg.get_p2().get_xy()
    if not math.isclose(ptl[1], 0.0, rel_tol=0.0,
                        abs_tol=TOL_ACC) or not math.isclose(
            ptr[1], 0.0, rel_tol=0.0, abs_tol=TOL_ACC):
        # initial segment should be on x-axis
        print("problem with initial hull segment")
        raise ValueError
    if ptl[0] > ptr[0]:
        tmp_arr = ptl
        ptl = ptr
        ptr = tmp_arr
    p_left = Point(*ptl)
    p_right = Point(*ptr)
    l_hull.append([p_left, None])
    r_hull.append([p_right, None])
    return l_hull,r_hull


def hull_seg_start(l_hull,r_hull,seg):
    # this helps start the convex hull on the respective first upward-facing
    # side of the starting triangle; this will be called ~right after calling
    # hull_init();
    # this changes mutables l_hull and r_hull in place;
    # there can be a special case, because the "unfolding" may not pivot off
    # the base segment endpoint if the base subinterval (of base triangle) is
    # a proper subinterval of the base triangle's base side (on x axis)--this
    # requires ~kludging a segment from end of subinterval to either end of
    # base triangle's base side
    p_left = l_hull[0][0]   # expects l_hull just one point
    ptl = p_left.get_x()
    p_right = r_hull[0][0]  # expects r_hull just one point
    ptr = p_right.get_x()
    p1, p2 = seg.get_p1p2() # base triangle side endpoints;
    # it's possible the base segment goes all the way to triangle edge--if not,
    # requires special initialization
    if p1==p_left:  # unfolding off l_hull point
        r_hull[0][1] = get_angle(p_right,p2)
        r_hull.append([p2,None])
    elif p1==p_right:   # unfolding off r_hull point
        l_hull[0][1] = get_angle(p_left,p2)
        l_hull.append([p2,None])
    elif p2==p_left:    # unfolding off l_hull point
        r_hull[0][1] = get_angle(p_right,p1)
        r_hull.append([p1,None])
    elif p2==p_right:   # unfolding off r_hull point
        l_hull[0][1] = get_angle(p_left,p1)
        l_hull.append([p1,None])
    else:
        if math.isclose(p1.get_y(), 0.0, rel_tol=0.0, abs_tol=TOL_ACC):
            if p1.get_x() < ptl:
                l_hull[0][1] = np.pi
                l_hull.append([p1, None])
                r_hull[0][1] = get_angle(p_right, p2)
                r_hull.append([p2, None])
            elif p1.get_x() > ptr:
                r_hull[0][1] = 0.0
                r_hull.append([p1, None])
                l_hull[0][1] = get_angle(p_left, p2)
                l_hull.append([p2, None])
            else:
                print("problem with initial segments")
                raise ValueError
        elif math.isclose(p2.get_y(), 0.0, rel_tol=0.0, abs_tol=TOL_ACC):
            if p2.get_x() < ptl:
                l_hull[0][1] = np.pi
                l_hull.append([p2, None])
                r_hull[0][1] = get_angle(p_right, p1)
                r_hull.append([p1, None])
            elif p2.get_x() > ptr:
                r_hull[0][1] = 0.0
                r_hull.append([p2, None])
                l_hull[0][1] = get_angle(p_left, p1)
                l_hull.append([p1, None])
            else:
                print("problem with initial segments")
                raise ValueError
        else:
            print("problem with initial segments")
            raise ValueError


def get_angle(p1,p2):
    # determine the absolute angle of a line segment from p1 to p2 (p1
    # imagined at the origin);
    # returns an angle in (-pi,pi)
    return math.atan2(p2.get_y()-p1.get_y(),p2.get_x()-p1.get_x())


def angle_diff(base_ang,test_ang):
    # which side of an oriented line determined by base_ang is a vector at angle
    # test_ang on?; returns a value in [-pi,pi], w/ (0,pi) being on "CCW" side,
    # and (-pi,0) on "CW" side of a ray at angle base_ang
    raw_dff = test_ang-base_ang
    return np.mod(raw_dff+np.pi,2*np.pi)-np.pi


def update_hull_side(on_side,on_hull,nx_pt):
    # updates the given hull side and returns a new list as the modified hull
    # side (output hull side is a deepcopy)
    out_hull = copy.deepcopy(on_hull)
    for jj in range(len(out_hull) - 1):
        theta = get_angle(out_hull[jj][0], nx_pt)
        if angle_diff(theta, out_hull[jj][1]) <= 0:
            # if theta >= out_hull[jj][1]:
            if on_side == "left":
                pass
            else:  # on_side=="right"
                out_hull[jj][1] = theta
                out_hull = out_hull[0:jj + 1] + [[nx_pt, None]]
                break
        else:  # theta < out_hull[jj][1]
            if on_side == "left":
                out_hull[jj][1] = theta
                out_hull = out_hull[0:jj + 1] + [[nx_pt, None]]
                break
            else:
                pass
    else:
        out_hull[-1][1] = get_angle(out_hull[-1][0], nx_pt)
        out_hull.append([nx_pt, None])
    return out_hull


def seg_int_pts(p11,p12,p21,p22):
    # points input form of segment intersection test; segment 1 is p11,p12,
    # and segment 2 is p21,p22
    if p11==p21 or p11==p22 or p12==p21 or p12==p22:
        return True
    tot_pts = [p11.get_xy(),p12.get_xy(),p21.get_xy(),p22.get_xy()]
    mn_vs = [list(map(lambda a,b:b-a,tot_pts[0+2*ii],tot_pts[1+2*ii])) for ii in
            range(2)]   # these are the p11-to-p12 and p21-to-p22 segment
    # vectors
    nrms = [math.sqrt(tup[0]**2+tup[1]**2) for tup in mn_vs]
    sgns = [0.0,0.0]
    axs_tst = [[False,tuple()],[False,tuple()]]
    for ii in range(2): # each pass is just a dot product
        # dot product with perpendicular:
        tmp = list(map(lambda a,b: a*b,list(map(lambda a,b:b-a,tot_pts[0],
                tot_pts[ii+2])),[mn_vs[0][1]/nrms[0],-mn_vs[0][0]/nrms[0]]))
        sgns[ii] = tmp[0]+tmp[1]
        if abs(sgns[ii])<TOL_ACC:   # point is effectively on axis of mn_vs[0]:
            # projection onto unit perpendicular gives distance away from
            # segment; if distance is within the similar "tube" as in
            # TriangleSegment's in_interior() method, then consider the point
            # on axis
            axs_tst[ii][0] = True
            # dot product with parallel:
            tmp_int = list(map(lambda a,b: a*b,list(map(lambda a,b:b-a,
                tot_pts[0],tot_pts[ii+2])),[mn_vs[0][0]/nrms[0],
                                            mn_vs[0][1]/nrms[0]]))
            tmp_pos = tmp_int[0]+tmp_int[1]
            axs_tst[ii][1] = (tmp_pos<-TOL_ACC,tmp_pos>nrms[0]+TOL_ACC)
            if not axs_tst[ii][1][0] and not axs_tst[ii][1][1]:
                return True # the point is effectively interior to the segment
    if axs_tst[0][0] or axs_tst[1][0]:
        if axs_tst[0][0] != axs_tst[1][0]:
            return False    # one is on the axis, and not an interior point,
            # and the other is not on the axis--no way to intesect
        else:   # both other endpoints are
            # on axis; check if they straddle the segment
            return not((axs_tst[0][1][0] and axs_tst[1][1][0]) or (axs_tst[
                0][1][1] and axs_tst[1][1][1]))
    elif (sgns[0]>0.0 and sgns[1]>0.0) or (sgns[0]<0.0 and sgns[1]<0.0):
        return False
    axs_tst = [[False, tuple()], [False, tuple()]]
    for ii in range(2):
        tmp = list(map(lambda a,b: a*b,list(map(lambda a,b:b-a,tot_pts[2],
                tot_pts[ii])),[mn_vs[1][1]/nrms[1],-mn_vs[1][0]/nrms[1]]))
        sgns[ii] = tmp[0]+tmp[1]
        if abs(sgns[ii]) < TOL_ACC:  # same as comment above
            axs_tst[ii][0] = True
            tmp_int = list(map(lambda a, b: a * b, list(map(lambda a,b: b-a,
                    tot_pts[2],tot_pts[ii])), [mn_vs[1][0]/nrms[1],
                                               mn_vs[1][1]/nrms[1]]))
            tmp_pos = tmp_int[0] + tmp_int[1]
            axs_tst[ii][1] = (tmp_pos < -TOL_ACC, tmp_pos > nrms[1] + TOL_ACC)
            if not axs_tst[ii][1][0] and not axs_tst[ii][1][1]:
                return True  # the point is effectively interior to the segment
    if axs_tst[0][0] or axs_tst[1][0]:
        if axs_tst[0][0] != axs_tst[1][0]:
            return False  # one is on the axis, and not an interior point,
            # and the other is not on the axis--no way to intesect
            # (case of all on mutally same axis has been done above)
    elif (sgns[0] > 0.0 and sgns[1] > 0.0) or (sgns[0] < 0.0 and sgns[1] < 0.0):
        return False
    return True


def on_side_check(l_hull,r_hull,nxt_seg):
    # checks which convex hull side the next segment unfolds from; returns
    # None if no unfolding appears to have occurred (points don't intersect)
    p1,p2 = nxt_seg.get_p1p2()
    if p1==l_hull[-1][0]: # unfolding hinge on left; modifying right
        on_side = "right"
        nx_pt = p2
    elif p2==l_hull[-1][0]: # unfolding hinge on left; modifying right
        on_side = "right"
        nx_pt = p1
    elif p1==r_hull[-1][0]: # unfolding hinge on right; modifying left
        on_side = "left"
        nx_pt = p2
    elif p2==r_hull[-1][0]: # unfolding hinge on right; modifying left
        on_side = "left"
        nx_pt = p1
    else:
        on_side = None
        nx_pt = None
    return (on_side,nx_pt)


def check_hull_intersect(l_hull,r_hull,on_side):
    # main function for handling intersection testing between the l and r
    # convex "hallway" hulls; this assumes only the last segment in the on_side
    # hull ("left" or "right") needs to be checked against all segments in
    # the other hull (for intersection);
    # note, hull lists have most recent segments first (0 index);
    # returns True if an intersection is found; False otherwise
    if on_side=="right":
        new_hull = r_hull
        old_hull = l_hull
    else:
        new_hull = l_hull
        old_hull = r_hull
    ns_vc = make_vec(new_hull[-1][0],new_hull[-2][0],norm=False)
    ns_nm = np.linalg.norm(ns_vc)
    un_vc = ns_vc / ns_nm
    hull_pt_ls = [tup[0].get_xy() for tup in old_hull]  # list of hull points
    # (lists) [[x1,y1],[x2,y2],...]
    bs_pt_bc = [new_hull[-1][0].get_xy()] * len(hull_pt_ls)  # broadcast new
    # segment "base point"
    chk_arr = np.array(hull_pt_ls)-np.array(bs_pt_bc)  # prep for broadcast dot
    # product
    prj_tot = np.sum(chk_arr * un_vc, axis=1)  # projection of each hull
    # point onto new segment line, origin at new_seg_p1
    lim_arr = [(prj_tot[ii] < 0.0 - TOL_ACC, prj_tot[ii] > ns_nm + TOL_ACC) for
               ii in range(len(prj_tot))]
    pnt_prs = [(ii,ii+1) for ii in range(len(prj_tot)-1) if not
        (lim_arr[ii][0] and lim_arr[ii+1][0]) and not (lim_arr[ii][1]
                                                   and lim_arr[ii+1][1])]
    for pr in pnt_prs:
        int_fnd = seg_int_pts(old_hull[pr[0]][0],old_hull[pr[1]][0],
                    new_hull[-1][0],new_hull[-2][0])
        if int_fnd:
            break
    else:
        return False
    return True


def plot_sides(side_list,clr='red',mrk=False):
    # for debugging; plots a list of triangle side objects
    for side in side_list:
        plt_lin(*side.get_p1p2(),clr=clr,mrk=mrk)


def plot_hulls(l_hull,r_hull,clr='green'):
    # for debugging; can plot convex hulls
    sl_left = []
    for jj in range(len(l_hull) - 1):
        sl_left.append(TriangleSide(l_hull[jj][0], l_hull[jj + 1][0]))
    sl_right = []
    for jj in range(len(r_hull) - 1):
        sl_right.append(TriangleSide(r_hull[jj][0], r_hull[jj + 1][0]))
    plot_sides(sl_left, clr=clr, mrk=True)
    plot_sides(sl_right, clr=clr, mrk=True)




###
# MAIN start
###

if __name__=="__main__":

    ### enter a pickle filename (.pkl), for saving vertex arcs data to
    out_pkl_fle = ""

    ### set the base side length for the base triangle (base side lies on x-axis,
    # with left end at origin)
    largest_side_length = 2.0   # a standard side length; this is not critical

    ### set the angles for the triangle
    angs = np.array([2/10,1/10,7/10])   # eg for 1,2,7 triangle

    ### set the interval on the base side that contains relevant trajectory foot
    # points for determining boundaries of discontinuity; set to [a,b],
    # 0<=a<b<=largest_side_length; defaults to [0,largest_side_length] if set
    # to None; note, sub_bas_int should always have lower x point first
    sub_bas_int = [0.25,0.55]

    ### set the radial interval that relevant trajectories must lie in for
    # determining boundaries of discontinuity (this is used in both the
    # unfolding process (discarding any side that lies outside a spread cone--ie a
    # ray w/ lower angle emanating from right end of sub_bas_int, and ray w/
    # higher angle emanating from left end of sub_bas_int), and in vetting the
    # final list of vertices and corresponding base intervals)
    radial_target = RadialInterval(2.1798,2.2417)

    ### set depth of unfoldings
    num_unfold = 10

    start_time = time.time()

    base_tri = triangle_from_angles(*(np.pi*angs),largest_side_length)

    if sub_bas_int is None:
        base_subinterval = base_tri.get_side(2)
    else:
        if sub_bas_int[0]>sub_bas_int[1]:
            sub_bas_int.reverse()
        base_subinterval = TriangleSide(Point(sub_bas_int[0],0.0),
                                        Point(sub_bas_int[1],0.0))

    edges_lst = []  # this will be a list of lists; the j'th entry corresponds
    # to edges derived from the j-1'th unfolding; j=0 corresponds to a single
    # edge, the base triangle's base; j=1 corresponds to base triangle's two
    # non-base sides; j=2 corresponds to all 4 sides of both unfoldings (
    # provided the sides are ~co-linear through their parents to the base
    # side); the j'th entry is itself a list, of tuples, w/ tuple format,
    # (index of parent edge in list of tuples at edges_lst[j-1], triangle
    # side object)

    hulls_lst = []  # this will be a list of tuples, [left_hull,right_hull],
    # corresponding to segment paths that begin w/ respectively indexed most
    # recent (last) segment sublist in edges_lst--eg edges_lst[-1][0] will
    # correspond to hulls_lst[0]'s pair of convex hulls, etc.
    # l_hull and r_hull each will contain sublists of type [point,angle],
    # where point is the start point of that facet, and angle is the angle of
    # vector from start point to the next facet's start point (angle=None if
    # none)


    '''
    notes:
    
    edges_lst[0]: [(-1,base side of base_tri)]
    edges_lst[1]: [(0,side 0 of base tri), (0, side 1 of base tri)]	# 1st
    indices refer to index of "parent" side in the list one index back in
    edges_lst
    
    lead_tri_lst: [(0,1),base_tri]  # the pair in the tuple are the indices in  
    edges_lst[1] corresponding (in order) to the "outer-facing" edges in the 
    triangle to be unfolded (base_tri)
    
    to unfold a triangle, we'll have [(ix_0,ix_1),tri]
        new_tri = tri.unfold(0) # unfold along side indexed 0
        add new_tri side 0 to edges_lst, ie tuple (ix_0,edge 0), and record
            index, ix_a
        add new_tri side 1 to edges_lst, ie tuple (ix_0,edge 1), and record
            index, ix_b
        add new_tri itself to tri_lst, ie tuple ((ix_a, ix_b), new_tri)
        <repeat for unfold(1)>
    '''

    # plot initialization; do here if wanting debugging plot capability
    # in main loop
    #fig = plt.figure()
    #ax = fig.add_subplot(1, 1, 1)

    #*****
    # initialize main loop lists
    #*****

    sid_nms = base_tri.get_side_nums()
    edges_lst.append([(-1,sid_nms[2],base_subinterval)])   # first entry,
    # just base side of base triangle or subinterval of it
    edges_lst.append([(0,sid_nms[0],base_tri.get_side(0)),(0,
                                sid_nms[1],base_tri.get_side(1))])  # second
    # entry, both upward-facing sides of initial triangle

    for tup in edges_lst[1]:
        l_hull,r_hull = hull_init(base_subinterval)
        hull_seg_start(l_hull, r_hull, tup[2])
        hulls_lst.append([l_hull,r_hull])

    lead_tri_lst = [((0,1),base_tri)]

    vrt_vis_bas_lst = []    # list of 2-tuples: (vertex point,visible interval
    # on base triangle's base side (along x-axis))--each element in the (2-)
    # tuple is a 2-element list of floats; note, the visible interval should
    # correctly be both within the base_subinterval's segment, and all angles
    # from the visible interval to the vertex should be within radial_target's
    # limits

    # check base triangle vertex "visibility" to base side
    cur_vrt = base_tri.get_vertex(1)
    vis_arc = get_vertex_vis(cur_vrt,[base_subinterval])
    tmp_arc = vtx_arc_int_con(cur_vrt, vis_arc[0],base_subinterval,
                    radial_target)
    if tmp_arc is not None:
        vrt_vis_bas_lst.append(([cur_vrt.get_x(), cur_vrt.get_y()],
                    tmp_arc))


    #*****
    # primary unfolding iteration loop
    #*****
    for jj in range(num_unfold): # number of unfolding levels to commit

        tmp_tri_lst = []
        new_edges = []
        new_hulls = []
        for tup in lead_tri_lst:   # for each "leading edge" triangle
            triangle = tup[1]
            tp_ix = tup[0]
            for side_ix in range(2):    # for each of 2 sides of the triangle
                # to unfold
                if tp_ix[side_ix]==-1:  # possible that a side is not visible;
                    # skip
                    continue
                sde_lst = fetch_edges(tp_ix[side_ix], edges_lst)    # sde_lst
                # will be a list of edges in specific order--first element is
                # latest edge unfolded, while last element is base edge of base
                # triangle (or subinterval)
                new_tri = triangle.unfold(side_ix)  # unfold along side of
                # side_ix
                
                #*****
                # check vertex visibility
                #*****
                cur_vrt = new_tri.get_vertex(1)    # vertex to check visibility of
                # note, vertices can appear multiple times in this loop, but the
                # result of different unfoldings; each unfolding will present a
                # unique path of past segments to check vertex visibility through
                # visibility from base of base triangle
                vis_arc = get_vertex_vis(cur_vrt,sde_lst)
                if vis_arc is None: # this can be None--a vertex
                    # could be created from a "legitimate" side's unfolding,
                    # but that vertex could be ~out of sight through previous
                    # unfoldings to the base; otoh, it could be an error
                    pass
                elif len(vis_arc)==1:
                    cur_arc = vis_arc[0]
                    tmp_als = cur_arc.get_angs()
                    if tmp_als[0]>tmp_als[1]:
                        tmp_als[1] += np.pi
                    if abs(tmp_als[1]-tmp_als[0]) < TOL_ACC:
                        # check radial interval cur_arc is not trivial--it's
                        # possible the vertex is eg on the x-axis, which can produce
                        # a ~false, trivial radial interval for the sightline
                        pass
                    elif tmp_als[0] <= np.pi or tmp_als[1] >2*np.pi:
                        # visible arc (looking down from vertex
                        # "above") should be in [pi,2*pi]
                        print("problem with vertex visible arc")
                        raise ValueError
                    else:
                        # convert visible arc to interval on base triangle's
                        # base and store in vrt_vis_bas_lst list; after paring
                        # through the radial_target, the arc can be None,
                        # so don't save it if that's the case
                        tmp_arc = vtx_arc_int_con(cur_vrt, vis_arc[0],
                                        base_subinterval, radial_target)
                        if tmp_arc is not None:
                            vrt_vis_bas_lst.append(([cur_vrt.get_x(),
                                            cur_vrt.get_y()], tmp_arc))
                else:   # ie vis_arc contains more than one arc
                    # if multiple arcs are visible from cur_vrt, then
                    # something's probably wrong
                    print("problem with determining vertex visible arc")
                    # DEBUG:
                    plot_sides(sde_lst)
                    plt.scatter(*cur_vrt.get_np_arr(),marker='*')
                    plt.show()
                    # END DEBUG
                    raise ValueError
                
                #*****
                # check side visibility
                #*****
                side_nms = new_tri.get_side_nums()
                ne_ix = [-1,-1]
                for ii in range(2): # check each new side of newly unfolded
                    # triangle
                    tmp_sde = new_tri.get_side(ii)
                    # is the new edge visible?
                    # is it within radial target "cone" (if applicable)?
                    if radial_target is not None:
                        angs = radial_target.get_angs()
                        pt_ls = base_subinterval.get_p1p2_int()
                        bas_int = [pt_ls[0].get_x(),pt_ls[1].get_x()]
                        if bas_int[0]>bas_int[1]:
                            bas_int.reverse()
                        # check if new side is completely on one side or the
                        # other of the inclusion cone
                        on_rt = side_check_sor(tmp_sde,angs[0],bas_int[1],'right')
                        on_lt = side_check_sor(tmp_sde,angs[1],bas_int[0],'left')
                        if on_rt or on_lt:
                            continue
                    # update hulls
                    l_hull,r_hull = hulls_lst[tp_ix[side_ix]]
                    on_side, nx_pt = on_side_check(l_hull, r_hull, tmp_sde)
                    if on_side is None:  # does not appear to be a legitimate
                        # unfolding
                        print("problem with unfolding")
                        raise ValueError
                    if on_side == "left":
                        l_hull = update_hull_side(on_side, l_hull, nx_pt)
                        if len(l_hull) > 2: # check for spiraling
                            if angle_diff(l_hull[0][1], l_hull[-2][1])<0.0:
                                print("problem with convex hull spiraling")
                                # DEBUG:
                                plot_sides([tmp_sde]+sde_lst,clr='red')
                                plot_hulls(l_hull,r_hull,clr='green')
                                ax.set_aspect('equal')
                                plt.show()
                                tt_ls = [tmp_sde]+sde_lst
                                # END DEBUG
                                raise ValueError
                    else:
                        r_hull = update_hull_side(on_side, r_hull, nx_pt)
                        if len(r_hull) > 2: # check for spiraling
                            if angle_diff(r_hull[0][1], r_hull[-2][1])>0.0:
                                print("problem with convex hull spiraling")
                                raise ValueError

                    # check corner case of a single hull side "flipping" to
                    # below x axis (can do this, and still have hulls not
                    # intersect)
                    flip_hull = False
                    if on_side=="left":
                        nw_sg = l_hull[-1]  # latest pt in hull pts
                        if (nw_sg[0].get_y()<0.0 and
                            nw_sg[0].get_x()>base_subinterval.get_p1().get_x()):
                            flip_hull = True
                    else:
                        nw_sg = r_hull[-1]  # latest pt in hull pts
                        if (nw_sg[0].get_y()<0.0 and
                            nw_sg[0].get_x()<base_subinterval.get_p2().get_x()):
                            flip_hull = True
                    if not flip_hull:
                        # check if hulls intersect
                        chi = check_hull_intersect(l_hull,r_hull,on_side) #
                        # check for intersection of l and r convex hulls;
                        # this will be True if no line exists between the
                        # unfolded paths, and False otherwise
                        if not chi:
                            ne_ix[ii] = len(new_edges)
                            new_edges.append((tup[0][side_ix],side_nms[ii],
                                              tmp_sde))
                            new_hulls.append([l_hull,r_hull])
                        else:
                            pass
                if ne_ix != [-1,-1]:    # don't include if both edges are not
                    # visible (can both edges ever *not* be visible?)
                    tmp_tri_lst.append((tuple(ne_ix),new_tri))
        print("triangles leading edge no. %s: %s total triangles" % (jj,
                                                           len(lead_tri_lst)))
        print("number of edges in edges_lst: %s" % sum([len(x) for x in
              edges_lst]))

        edges_lst.append(new_edges)
        hulls_lst = new_hulls
        lead_tri_lst = tmp_tri_lst
        # note, at this point, lead_tri_lst has tuples ((ix_0,ix_1),tri_obj),
        # where the ix_'s reference indices in the list at the end (just added) of
        # edges_lst

        #for tup in tmp_tri_lst: # DEBUG
            # for plotting / checking new triangles at this iteration
            #tup[1].plot_triangle()

    print("number of arcs found: %s" % len(vrt_vis_bas_lst))
    print("total time: %s" % (time.time() - start_time))

    ### to save vertex boundary arcs data (for possible use with other files):
    if out_pkl_fle != "":
        with open(out_pkl_fle,'wb') as fp:
            pickle.dump(vrt_vis_bas_lst,fp)


    ### optional plotting and visualization section:

    #****
    # display base triangle and points of unfolded vertices visible through
    # the (x,theta) window
    #****

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    # highlight, with markers, which vertices in the total unfolding are visible
    # from the base triangle base side's sub_bas_int range, and the radial_target
    # span of acceptable trajectory angles
    # angles
    for tup in vrt_vis_bas_lst:
        if tup[1] is not None:
            plt.scatter(*tup[0],marker='o',s=400,color='red',facecolors='none')

    # highlight the base triangle on the plot
    base_tri.plot_triangle(clr='red')

    ax.set_aspect('equal')
    plt.show()


    #*****
    # display discontinuity boundaries / boundary arcs
    #*****

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    arc_res = 500    # resolution level for plotting the arcs
    tmp_p1,tmp_p2 = base_subinterval.get_p1p2()
    x_msh = abs(tmp_p2.get_x()-tmp_p1.get_x())/arc_res # x mesh size
    if radial_target is not None:
        tmp_a1,tmp_a2 = radial_target.get_angs()
    else:
        tmp_a1,tmp_a2 = 0.0,np.pi
    if tmp_a2 < tmp_a1:
        tmp_a2 += 2*np.pi
    y_msh = (tmp_a2-tmp_a1)/arc_res # angular mesh size (y coordinate in arc
    # plots)

    for tup in vrt_vis_bas_lst:
        if tup[1] is None:
            continue
        pt_ls = []  # for holding all point tuples (x,y) on the arc
        vt_pt = tup[0]  # x,y of vertex point
        bs_iv = tup[1]  # [a,b] of base_subinterval subsegment vertex is visible
        # from through all the unfoldings
        ang_l = np.arctan2(vt_pt[1]-0,vt_pt[0]-bs_iv[0])    # angle from
        # leftmost point of bs_iv
        ang_r = np.arctan2(vt_pt[1]-0,vt_pt[0]-bs_iv[1])    # angle from
        # rightmost point of bs_iv
        phs_vec = np.array([bs_iv[1]-bs_iv[0],ang_r-ang_l]) # vector in x,
        # ang phase space between the two interval-angle limits
        ell_cns = (phs_vec[0]/x_msh)**2+(phs_vec[1]/y_msh)**2 # for estimating
        # the number of mesh points along phs_vec, given aspect ratio re x_msh
        # and y_msh
        x_inc = np.sqrt(phs_vec[0]**2 / ell_cns)  # x mesh after accounting for
        # aspect ratio
        x_vals = np.arange(bs_iv[0],bs_iv[1],x_inc)
        a_vals = [np.arctan2(vt_pt[1],vt_pt[0]-x) for x in x_vals]
        plt.plot(x_vals,a_vals,color='blue')

    plt.gca().invert_yaxis()
    ax.set_aspect('equal')
    plt.show()


