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

load "data/AF-P11362-F1_clean.pdb", protein
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

load "data/AF-P11362-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [3539,2961,3110,3487,3988,3989,3492,2903,2905,2910,2959,3509,3513,3530,2904,2907,2908,3126,3127,4069,2962,3960,4064,4066,4067,3370] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.329,0.278,0.702]
select surf_pocket2, protein and id [1857,1856,1858,4562,2334,2340,2341,1866,1869,1870,1872,1874,1875,1876,1877,1878,1853,4221,4223,4561,4218,4560,2101,2102,4105,4525,4527,4528,4523,4526,4493,4495,4521,2343,4524,4211,4244,4497,4496,4529] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.698,0.361,0.902]
select surf_pocket3, protein and id [2926,2923,3126,3127,4071,4072,4073,2962,3257,3125,3233,3122,3124,2946,2934,2943,2944,2945,3132,2951,4076,4078,3221,3222,3223,3228,3230,2929,2935,2939,2941,3143] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.639]
select surf_pocket4, protein and id [4491,4156,4157,4158,4085,4086,4087,4088,4089,4172,966,1357,4475,4467,1356,1355,4440,4447,4438,1370,1227,1229,1230] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.545]
select surf_pocket5, protein and id [2321,2307,2375,2377,2378,2379,2380,2381,1916,2308,1935,1936,2288,2397,2036,2055] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.353,0.278]
select surf_pocket6, protein and id [4024,4025,4026,4027,3744,3747,3769,3719,3580,3701,3745,3767,3800,3775] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.729,0.361]
select surf_pocket7, protein and id [3870,3869,3842,3841,3868,4825,4828,4821,4824,4275,4277,4290,4823] 
set surface_color,  pcol7, surf_pocket7 
   

deselect

orient
