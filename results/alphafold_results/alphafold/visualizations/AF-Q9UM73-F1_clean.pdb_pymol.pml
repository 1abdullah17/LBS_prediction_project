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

load "data/AF-Q9UM73-F1_clean.pdb", protein
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

load "data/AF-Q9UM73-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5142,5155,5157,5959,5961,4411,8830,8831,4410,8790,4436,4437,5956,5957,5102,4416,4434,4435,4438,4424,4428,4429,4407,4409,5984,8833,5170,5168,5166,8726,8730,8920,8757,8759,8760,8761,8725,8753] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.333,0.702]
select surf_pocket2, protein and id [6953,6956,7010,7141,7997,7998,7547,7544,7548,7531,7568,7572,7561,7562,7565,6957,8084,8085,8089,8093,7970,7525,7008,7011,6952,6954,6958,6959,6960,6973,7157,7158,6966,7156] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.443,0.361,0.902]
select surf_pocket3, protein and id [5965,5967,4453,4455,5968,4470,4466,8751,8626,8629,8614,8619,8624,8622,8623,5975,5976,4329,4330,4332,6011,5964,5970,4314,4316,5966,8644,8651,8652,5951,5931,5932,8716,8719,8750,8752] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.463,0.278,0.702]
select surf_pocket4, protein and id [3819,3525,3527,3530,2944,2948,3536,3538,3545,3211,3213,3529,3531,3534,3812,3528,3532,3524] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.741,0.361,0.902]
select surf_pocket5, protein and id [6976,6979,6991,8107,8108,8106,8115,8116,8114,8140,7927,7947,7983,7984,8089,8090,8091,8092,8093,7943,8109,8096,8105] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.694,0.278,0.702]
select surf_pocket6, protein and id [150,3874,3888,3889,3890,147,3878,2683,2684,3856,2661,2662,3876,3877,3116,2655] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.757]
select surf_pocket7, protein and id [2731,2734,2698,2693,2710,2738,2689,2755,2757,2758,2668,2629,2633,2691,2635,2638,2660,2664] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.475]
select surf_pocket8, protein and id [109,148,150,3888,3889,3890,142,144,146,1121,145,147,181,182,1120,2661,111,117,119] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.459]
select surf_pocket9, protein and id [5504,5406,5435,5437,5486,5503,5403,5501,3995,3996,5483,5490,5485] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.318,0.278]
select surf_pocket10, protein and id [1797,1819,1825,1827,1831,1836,1780,1787,1791,1793,1837,1840,1801,2284,2275,1859,1861,2219,2221] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.561,0.361]
select surf_pocket11, protein and id [3321,3323,3324,3700,3701,3686,3713,3717,3719,3426,3307,3309,3296,3381] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.553,0.278]
select surf_pocket12, protein and id [2515,2529,2530,3819,3548,3550,3552,3554,3213,3558,2531,2514] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.859,0.361]
select surf_pocket13, protein and id [137,138,172,173,174,3121,3122,3127,108,162] 
set surface_color,  pcol13, surf_pocket13 
   

deselect

orient
