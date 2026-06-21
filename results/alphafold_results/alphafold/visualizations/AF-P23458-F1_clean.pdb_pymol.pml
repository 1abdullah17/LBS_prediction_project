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

load "data/AF-P23458-F1_clean.pdb", protein
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

load "data/AF-P23458-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5154,5156,6400,6403,4199,5155,4198,4818,5151,4801,4802,5152,5180,4791,4793,4369,4370,4371,4133,4136,4770,4776,4771,4783,5463,5465,5468,5469,5733,5678,5431,5423,5425,5815,5816,6414,6399,6396,6398,6397,4825,6607,4134,4139,6430,6609,6431,4137,6624,6625,5732,5734,6380,6383,6378,5746,5751,5754,5771,4822,4866,5163,5307,4388,5312,5314,5318,5320,4751,5164,4387,5179,5121,5122,5127,5140,5459,5126,5435,5436,5437,4146,4154,4155,4156,4157,4181,4158,4163,5403,5429,6995,5424,5426,6389,4142,6387,6993,6994,6390,4162,5315,5323,5324,5327,5440,4386] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.392,0.702]
select surf_pocket2, protein and id [2899,2923,2875,2884,2885,2887,2891,2897,217,2873,2890,3175,3180,3506,3508,3462,2979,2980,319,320,343,3217,3181,3171,3177,3215,3213,3226,2910,2912,2914,2915,3583,2729,2756,2757,2727,2769,2771,350,351,353,348,352,2759,2760,2724,2763,2796,218,482,2880,480,204,3483,215,210,241,346,247,220,222,224,223,248,249,255,662,484,485,2700,466,661,367,375,377,379,2701,2702,2795,3882,2971,3883,2793] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.431,0.902]
select surf_pocket3, protein and id [4988,4989,2082,2083,2084,4964,4983,4985,4990,4991,4995,4997,6188,6220,6214,2432,2530,2521,2524,2525,2553,2555,2572,2554,2556,6153,2441,2442,2576,2575,2507,2509,2577,2457,2454,2456,6222,5017,2559,4961,4963,4965,4966,5253,4967] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.278,0.702]
select surf_pocket4, protein and id [7077,7081,7099,7102,7497,7133,6476,6477,7073,7075,6461,6472,6473,6474,6665,7056,7061,6922,6923,6532,6534,6535,7495,7496,7573,7574,7576,7055,7131,7132,7468,7577,7482,7579,7581,7583,7476,6480,6482,6479,7477,7479] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.431,0.361,0.902]
select surf_pocket5, protein and id [1977,1988,1976,1982,1947,1956,1963,1974,1975,1969,1971,1973,437,1059,1060,1061,1064,1076,395,396,398,435,629,632,1063,1024,434,458,460,645,1990,1991,637,992,1026,666,668,670] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.392,0.278,0.702]
select surf_pocket6, protein and id [6391,6392,6986,6987,7032,6992,4148,4149,4152,4190,4191,6435,6438,6709,6711,7031,6672,6686,6688,6544,6542,6543,6528,6693,6526,6418,6419,6421,6558,4165,4409,4410,4411,4424,4164,4156,5404,4169,5403,5405,6991,4433,4434,7010,4171,4162,5377,5406] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.576,0.361,0.902]
select surf_pocket7, protein and id [2801,2930,2931,2932,2825,2828,2901,2904,2905,2934,2938,2814,6105,6107,6108,6109,6110,2929,3810,6136,6164,6170,6171,6011,6015,6172,6100,6086,6093,6094,6096,6069,6075,6077] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.506,0.278,0.702]
select surf_pocket8, protein and id [79,1291,88,98,100,1415,1293,1295,1407,946,1282,1314,1315,94,947,983,1281,1287,1285,1312,1324,95,107,104,106,1447,1460,1462,1413,1463,1465,1007,1008,1010,1041,1446,1449,984,975,977,1006,1001] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.722,0.361,0.902]
select surf_pocket9, protein and id [3793,2943,2944,2950,3809,3660,3659,3603,3628,3629,3623,3141,3625,2946,3795,3601,2937,6057,6080,6081,3774,3784,3794] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.616,0.278,0.702]
select surf_pocket10, protein and id [7606,7765,7605,7814,7816,7437,7712,7714,7715,7713,7763,7778,7729,7731,7745,7760,8225,7739,7874,8226,8223,7893,7894,7848,7842,7851] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.863,0.361,0.902]
select surf_pocket11, protein and id [3871,2994,3873,3874,2992,330,3191,3199,882,888,883,3192,3008,3009,3010,838] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.671]
select surf_pocket12, protein and id [4849,4851,4852,4853,4854,4855,5724,5725,6331,5726,6604,5732,5734,5756,5758,6382,6383,6413,6414,4826,6594,6330,6593,5717,4818,6403] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.792]
select surf_pocket13, protein and id [1181,1214,1500,1503,1186,1492,1493,1496,1497,1498,1525,1527,1089,1077,1080,1082,1083,1191,1182,1109,1156,1110] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.561]
select surf_pocket14, protein and id [5089,3823,3824,3825,5566,5569,5573,5563,5565,5576,5571,5572,6106,3832,3833,6108,5598,6062,6086,5064,5599] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.647]
select surf_pocket15, protein and id [8563,7940,8547,8551,8556,7910,7913,7915,7916,7918,7914,7346,8553,8536,8545] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.447]
select surf_pocket16, protein and id [4189,4190,4191,4208,4211,6438,6439,6442,6543,4187,4188,4418,4394,4417,4185,4395,4105,4209,4106,4710] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.506]
select surf_pocket17, protein and id [692,766,767,768,742,1914,1915,1921,1922,1926,689,691,1920,1953,1962,1916] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.337]
select surf_pocket18, protein and id [3973,3954,2776,3932,3933,3971,2780,3859,3890,825,2783,3837,3838,3956,2407,3861,3967] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.361]
select surf_pocket19, protein and id [1773,1780,1306,1308,1302,1305,1339,1348,1344,1738,1742,1299,1352] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.337,0.278]
select surf_pocket20, protein and id [1665,1667,1675,1677,1421,1643,1705,1674,1671,1412,1363,1364,1379,1382,1417,1359] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.506,0.361]
select surf_pocket21, protein and id [3694,3518,3522,3540,3542,3557,3618,3640,3643,3651] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.447,0.278]
select surf_pocket22, protein and id [6501,6502,6503,6504,7474,7767,7769,7482,7583,7445,7471,7124,7470,6486,6492,6493,6495,7804] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.647,0.361]
select surf_pocket23, protein and id [1870,1733,1695,1697,1253,1224] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.561,0.278]
select surf_pocket24, protein and id [2008,2010,1996,2658,1986,435,434,2678,447,458,1991] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.792,0.361]
select surf_pocket25, protein and id [844,846,831,868,860,3974,847,830,3941,3939,3940,801,803,3973,826,823,824] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.671,0.278]
select surf_pocket26, protein and id [5908,5909,5895,5907,5922,5993,5997,6028,5663,5914,5906] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.863,0.902,0.361]
select surf_pocket27, protein and id [8476,8475,8320,8130] 
set surface_color,  pcol27, surf_pocket27 
   

deselect

orient
