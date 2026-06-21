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

load "data/AF-Q6P5Z2-F1_clean.pdb", protein
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

load "data/AF-Q6P5Z2-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [4304,3405,3406,4306,4309,4311,4312,3283,4229,3284,3401,3403,3404,3684,3801,3802,3388,4189,4190,4214,2847,2850,4313,3231,3233,3240,3241,3242,3244,3267,2840,2843,2846,3247,3249,3250,3252,3253,3265,3266,2844,2855,3683,3807,5496,5498,3219,3221,3226,3229,2817,5467,2826,4202,2821,3852,4203,3874,5448,5457,5462,5464,3220,3223,3845,4230,3819,3823,4326,4327,4328,4336,3525,3561,3562,3528,3529,3531,3597,3264,3499,3257,3261,3262,3263,3462,3417,3419,2870,2871,2875,2876,2873,3487,3504,2861,2863] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.310,0.702]
select surf_pocket2, protein and id [593,595,597,598,599,600,601,1068,622,584,583,585,592,2606,2607,561,564,2613,2549,629,631,1036,1038,1039,1034,1037,2548,2550,2551,2552,2516,2561,2563,1040,1041,1108,1103,1104,1106,1109,2639,2625,2630,1102] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.498,0.361,0.902]
select surf_pocket3, protein and id [2929,2935,2937,2939,3090,2984,1328,3060,2943,2867,2900,2897,2901,2906,3087,3088,1418,1303,1304,1392,1393,1394,1355,1336,2923,1356,1359,2924,3067] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.525,0.278,0.702]
select surf_pocket4, protein and id [3052,1192,1195,1230,2104,2105,2106,2109,2107,2117,2341,1794,2111,2116,2118,2366,2506,1240,1556,2508,2354,2349,2352,2353,1202,1203] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.851,0.361,0.902]
select surf_pocket5, protein and id [2498,2502,2505,2660,2643,2644,1144,1176,2665,1145,1175,1611,1614,1581,1602,1608,1609,1641,1640,1201,1167,1165,2501,1206,2507,1578,1607,2645,1101,1102,1100,1108,1109] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.600]
select surf_pocket6, protein and id [2136,1460,1493,1494,1495,2095,1490,1492,2156,2397,2399,2400,1514,1516,2094,2370,2077,2078,2381,2474,1459] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.596]
select surf_pocket7, protein and id [4473,4471,2827,2802,2810,4768,4817,4761,2806,2805,4737,2829,4773,2329,4779,4784,2787,2795] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.325]
select surf_pocket8, protein and id [4351,2918,4350,3516,3523,3551,3517,3553,4355,4360,4349,4348,2882,3527] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.478,0.361]
select surf_pocket9, protein and id [5330,5239,5241,5331,5243,5249,3990,3994,5336,3981,3982,5248,5247,3959,3987,3988] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.510,0.278]
select surf_pocket10, protein and id [3014,3019,3002,3003,3004,3006,3020,1247,1249,3022,3017,4835] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.831,0.361]
select surf_pocket11, protein and id [2844,2855,3249,3250,3252,3253,2860,2870,2875,2903,2906,2873,3487,2902] 
set surface_color,  pcol11, surf_pocket11 
   

deselect

orient
