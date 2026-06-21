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

load "data/AF-P04626-F1_clean.pdb", protein
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

load "data/AF-P04626-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5736,6123,5724,5733,5739,5741,5765,6094,5801,5761,5763,6200,6108,6205,6207,6208,6209,5185,5170,5692,5338,5205,5339,6122,5161,5203,5206,5160,5705,5708,5352,5355,5356,5349,5340,5518,5686,5687,6201,5154,5157,7222,7223,5153,5155,7216,7218,7220,5709] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.333,0.702]
select surf_pocket2, protein and id [3595,25,26,30,33,35,3609,210,3386,3606,3612,21,20,22,3591,3593,3594,3602,36,270,2224,274,269,2225,2131,2156,2158,2123,2073,2133,2142,2145,3379,3383,3156,3385,3154,3159,3180,3387,3396] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.443,0.361,0.902]
select surf_pocket3, protein and id [1816,1818,1819,1820,1783,1839,1840,1778,1781,1854,1763,1769,1774,1749,2088,258,6,447,444,445,446,433,437,438,619,1838,2082,2081] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.463,0.278,0.702]
select surf_pocket4, protein and id [2958,2961,2462,2463,2468,2469,2531,2535,2537,2464,3192,3193,3188,3191,91,3199,3178,3172,2466,2538,2713,2712,2976,2978] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.741,0.361,0.902]
select surf_pocket5, protein and id [36,270,2224,268,275,276,667,285,668,3159,3203,3201,43,55,51,45] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.694,0.278,0.702]
select surf_pocket6, protein and id [2436,2667,2669,2671,2237,2268,2233,2236,2452,2231,2479,2945,2232,2664,2252,2487,2437,2455,2412,2273,2414,2419,2281] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.757]
select surf_pocket7, protein and id [7309,5716,5718,5719,7233,5704,7256,7259,7268,7235,5327,7344,5701,7326,7325,7328,7305] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.475]
select surf_pocket8, protein and id [1800,1217,1214,1799,1789,1790,1791,1801,1199,1200,1201,1788,987,988,990] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.459]
select surf_pocket9, protein and id [3492,3468,3495,3503,3504,3505,3699,3677,3671,3668,3665,3667,3669,3473,3650] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.318,0.278]
select surf_pocket10, protein and id [3339,3947,3971,3337,3335,3336,3325,3327,3897,3901,3908,3891,3967,3946,3913,3338] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.561,0.361]
select surf_pocket11, protein and id [2452,2231,2479,2667,2222,2288,2945,2232,2664] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.553,0.278]
select surf_pocket12, protein and id [170,3218,63,148,107,116,144,136,138,142,103,3198,3220,3641,140,110,112] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.859,0.361]
select surf_pocket13, protein and id [6895,6897,6899,6894,6896,6475,6615,6744,6303,6326,6705] 
set surface_color,  pcol13, surf_pocket13 
   

deselect

orient
