# Netgen .geo input
# 
# $Id: ellipsoid.geo,v 1.1 2007/05/14 12:58:58 zlb Exp $
#
# Note: After generating .mesh file, convert it to ALBERT format (using
# phgExportALBERT), and append the following lines to the end of the ALBERT
# file:
#
#	curved boundaries:
#	1
#	((x-.5)*4)^2 + ((y-.5)*5)^2 + ((z-.5)*5)^2 - 1;
#	.5 + (x-.5) / sqrt(((x-.5)*4)^2 + ((y-.5)*5)^2 + ((z-.5)*5)^2);
#	.5 + (y-.5) / sqrt(((x-.5)*4)^2 + ((y-.5)*5)^2 + ((z-.5)*5)^2);
#	.5 + (z-.5) / sqrt(((x-.5)*4)^2 + ((y-.5)*5)^2 + ((z-.5)*5)^2)

algebraic3d

solid ellipse = ellipsoid(0.5,0.5,0.5; 0.25,0,0; 0,0.2,0; 0,0,0.2) -bc=3;
solid cube = orthobrick(0, 0, 0; 1, 1, 1) and not ellipse -bc=1;

#-------------------------------------------------------------------------
# mesh (cube-ellipsoid) and ellipsoid
###tlo cube;
###tlo ellipse;
#-------------------------------------------------------------------------
# force aligning mesh faces with the plane (0.5,0.5,0.5;0,0,1)
solid cube0 = cube and plane(0.5,0.5,0.5; 0,0,1) -bc=0;
solid cube1 = cube and plane(0.5,0.5,0.5; 0,0,-1) -bc=0;
solid ellipse0 = ellipse and plane(0.5,0.5,0.5; 0,0,1) -bc=0;
solid ellipse1 = ellipse and plane(0.5,0.5,0.5; 0,0,-1) -bc=0;
tlo cube0;
tlo cube1;
tlo ellipse0;
tlo ellipse1;
#-------------------------------------------------------------------------
