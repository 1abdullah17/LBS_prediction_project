from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/AF-P53350-F1_clean.pdb", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta




# SAS points

load "data/AF-P53350-F1_clean.pdb_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [1209,1217,1316,1318,1320,1321,1322,1180,1343,1344,1335,1337,1345,798,1234,1231,1232,1233,1235,1236,797,799,1315,1336,1338,641,207,209,214,219,221,226,227,230,243,829,832,828,2656,1207,854,1179,1197,1476,838,861,862,863,890,366,264,378,379,380,775,535,536,395,241,512,472,500,433,503,471,794,781,566,533,640,248,234,237,239,245,240,208,211,262,760,509] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.298,0.702]
select surf_pocket2, protein and id [2990,2992,3467,3442,3444,3446,3443,2864,2868,3487,3002,2996,2999,3004,3007,3009,3011,3029,3400,3396,3399,3405,2288,2295,2697,2282,2285,2287,2276,2278,3445,2289,2266] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.533,0.361,0.902]
select surf_pocket3, protein and id [2704,2708,2710,2716,2720,3857,3694,3758,3759,3760,3856,3888,3887,3890,3894,2730,3373,2722,2723,2726,2728,2734,2738,2740,2744,3700,3738,3740,3742,3743,3891,2742,2746,2747,3713,3691,3739,3889,3365,3366,3910,3906,3908] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.565,0.278,0.702]
select surf_pocket4, protein and id [896,3804,3650,894,2316,895,2673,2679,898,3651,3652,2691,3649,3663,3790,2675] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.878]
select surf_pocket5, protein and id [1775,1490,1501,1503,1748,1749,853,878,881,1201,1204,1705,854,1194,1195,1198,1206,1475] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.533]
select surf_pocket6, protein and id [2769,2771,2774,2795,2797,2807,3254,3277,3278,3279,3280,3283,2809,3281,3321,3323,3325] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.490]
select surf_pocket7, protein and id [768,767,157,158,160,688,133,134,91,675,660,137,661,668,306,361] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.329,0.278]
select surf_pocket8, protein and id [1104,1595,1597,1599,2133,2125,1601,1099,1100,1107,1101,2117,1610,2113,2115,2119,1608] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.620,0.361]
select surf_pocket9, protein and id [2418,2381,2387,3804,895,2315,2408,3644,3647] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.631,0.278]
select surf_pocket10, protein and id [3300,3302,3339,3382,3337,3305,2874,2875,2845,2846,2872,2873,3003,3335] 
set surface_color,  pcol10, surf_pocket10 
   

deselect

orient
