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

load "data/AF-P43405-F1_clean.pdb", protein
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

load "data/AF-P43405-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [3449,3450,3529,3531,2652,2653,3537,3538,3539,3435,2895,3014,3015,3421,3423,3037,3039,3040,3042,2636,3020,3056,3065,3064,2511,2513,2514,2500,2465,2466,2468,2470,2472,2474,2475,2462,2459,2458,2460,2493] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.310,0.702]
select surf_pocket2, protein and id [1890,1891,1892,4467,1185,1218,1856,1858,929,930,931,152,153,154,1929,1932,868,1889,4466,4468,107,319,124,1189,1184,1867,1865,1866,1863,1896,1818,1820,1864,1801,1802,1900,4469,1191,1183,1190] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.498,0.361,0.902]
select surf_pocket3, protein and id [900,902,12,302,873,901,334,907,940,912,913,914,1,4,9,10,301,349,8,331,332,348] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.525,0.278,0.702]
select surf_pocket4, protein and id [4159,4161,4163,4190,4194,4401,4369,4370,4371,4372,4375,458,469,472,473,475,563,565,566,568,451,463,465,4406,4410,4446,4404,4162,4195] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.851,0.361,0.902]
select surf_pocket5, protein and id [1167,1169,1171,1160,4463,4459,3184,3185,3186,3191,3188,3223,3489,3490,3491,3197,2143,2160,2158,3194,1196,1198,4465,2142,2145,2147] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.600]
select surf_pocket6, protein and id [975,939,324,4436,331,936,972,937,335,1153,1154,4394,347,1126,1130,1132,1102,4390,1123,4438] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.596]
select surf_pocket7, protein and id [151,153,1939,1940,1941,188,185,1901,1954,1800,1801,1802] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.325]
select surf_pocket8, protein and id [58,74,76,266,785,786,628,648,647,651,652,654,658,649,688,627] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.478,0.361]
select surf_pocket9, protein and id [3633,3748,3750,3753,3773,3631,3634,3629,3630,3712,3715,3720,3718,3627,3651,3682,3684] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.510,0.278]
select surf_pocket10, protein and id [3548,3550,2894,2911,3014,3532,3533,3538,3561,2820] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.831,0.361]
select surf_pocket11, protein and id [956,2076,2079,921,1278,2065,2066,2073,2019,2031,2033,2021,2018,944,981,2029,2024,2028] 
set surface_color,  pcol11, surf_pocket11 
   

deselect

orient
