These are some basic geometric routines in R<sup>2</sup> for the unfolding of triangles in the context of mathematical billiards trajectories. See the associated [github.io page](link) for details.

In general, references to (x,theta) coordinates mean the following. A given triangle is oriented so that its longest side, its "base side," is along the x-axis, with the side's left endpoint at the origin. The other two sides of the triangle are in the +y halfplane. A billiard trajectory in this triangle is specified by its (x,theta) coordinates, where x is the trajectory's foot point on the triangle's base side, and theta is the angle made by the trajectory with respect to the +x axis, with theta in (0,pi).


## Unfoldings and discontinuity boundaries

The triangle_unfold_trajectory.py file will plot the unfolding of a given triangle along a given trajectory.

The triangle_discontinuity_arcs.py file will show discontinuity boundaries after a given number of unfoldings for a given triangle and a given phase space window in (x,theta) space. The triangle is specified by its interior angles, plus the length of its base side. The phase space window is specified by an x-axis subinterval within the base side segment, and a theta radial interval within (0,pi).

The triangle_discontinuity_arcs routine tracks all possible unfolding paths within the phase space window, and checks the sightline from vertex to base segment by using a convex hull reduction on the path of unfolded segments. This greatly speeds processing time.

![closeup of discontinuity boundary cells for unfolding of 3,4,5 triangle]([https://github.com/tmwine/mathematical_billiards/images/3_4_5_closeup.png](https://github.com/tmwine/mathematical_billiards/blob/main/images/3_4_5_closeup.png)?raw=true)

Though probably not as important, also included is the arcs_to_polygons.py file. This allows rendering the output of triangle_discontinuity_arcs as a set of adjacent polygons in (x,cotangent(theta)) phase space (after including the phase space boundary rectangle). The polygon data allows rendering the discontinuity cells as vector graphics files. The arcs_to_polygons routine requires the sweep_line module and its two AVL tree modules. These are available [nearby]([link](https://github.com/tmwine/proximal_BentleyOttmann_variation)).


## Bounce sequences and color plots

![closeup of 1,2,7 triangle, near angle pi/2]([https://github.com/tmwine/mathematical_billiards/images/1_2_7_closeup_4.png](https://github.com/tmwine/mathematical_billiards/blob/main/images/1_2_7_closeup_4.png)?raw=true)

To create the bounce sequence plots of the [github.io page](link), bounce sequences may first be generated with the billiard_bounce_sequences.jl file. This requires Julia, and its DynamicalBilliards package. The Julia routine runs from the command line as follows:
```
$ julia billiard_bounce_sequences.jl "[[x1,y1],[x2,y2],[x3,y3]]" low_theta:theta_step:high_theta low_x:x_step:high_x sequence_length
```
The triangle's coordinates are first specified (quotes required), where the triangle must be oriented as mentioned above: its longest side on the x-axis, its left endpoint, [x1,y1], required to be the origin, [0.0,0.0], and [x2,y2] must be [L,0.0], where L is the length of the longest side. The point [x3,y3] is the remaining vertex, y3>0. The next two ranges specify the phase space grid of all the initial trajectories. Finally, sequence_length is the desired length of the bounce sequences. The sequences are saved in a file you specify the filename of in the routine's header. One symbolic sequence is recorded per line. The Julia routine also displays the dimensions of the phase space grid--these are necessary for processing the symbolic sequence output file.

To process the output file of the Julia routine, three measures are available to convert symbolic sequences into real values:
- angular_sequence_measure.py--this is the simple, order-not-important measure of symbol proportionality; higher values correspond to more equal proportions
- lempel_ziv_measure.py--this is a variant of LZ77 compression; higher values correspond to greater compression; (this routine may be somewhat slow when fed a large number of sequences)
- MultAlpha_Fast_File_Batch.cpp--the matrix measure (needs compiling); it uses the companion file MiP_1_1_1_hk_4.bin, which has the basic 3-symbol 6x6 equiprobable matrices encoded in it (generated from [this](link) neighboring repo); to run, compile the .cpp file and assuming the executable is named MultAlpha_Fast_File_Batch, run from the command line as follows:
```
$ .<path to>/MultAlpha_Fast_File_Batch <path to>/MiP_1_1_1_hk_4.bin <path to>/measure_output_filename
```

All three measures will output measurement values to an intermediary text file, one value per line. To render the values as a color plot, you may run this short Octave script:
```
filename = "<path to>/measure_output_filename"
mx_vals = dlmread(filename,"\n");
mx_vals = reshape(mx_vals,num_pos,num_ang); % ie takes first num_pos elements and converts them to column 1, etc.; put another way, mx_vals now has "num_pos" rows, and "num_ang" columns
mx_vals = mx_vals';	% mx_vals now will have different positions along the columns, and different angles along the rows
figure;
imagesc(mx_vals);	% this does the color plot; Octave's default colormap is viridis
```
where num_ang and num_pos are the respective dimensions of the angle and x position ranges fed to the Julia bounce sequence script (dimensions which the Julia routine displays in its output).


## Dependencies

matplotlib, but only necessary if wanting intermediary visualizations
numpy

Julia is needed only if using the bounce sequence generating routine.

Octave is helpful for rendering bounce sequence value files as a color plot, but there are other applications that can do this.


## Licenses

This author's code is licensed under the MIT license.

The code in arcs_to_polygons uses AVL trees in its sweep_line routine from the neighboring [repo](https://github.com/tmwine/proximal_BentleyOttmann_variation). The AVL trees are from another 3rd party repo, under the MIT license. See the sweep_line repo for details.



