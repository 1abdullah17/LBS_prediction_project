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

load "data/AF-O60674-F1_clean.pdb", protein
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

load "data/AF-O60674-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5037,5038,5039,5635,6285,6286,6287,6288,6304,6305,6303,6506,6510,5589,5632,5633,5634,5063,5062,5209,5211,5047,4676,4678,5034,5035,5036,4673,6509,4646,4669,4670,6511,5339,5326,5327,5328,5332,5338,5368,5712,5329,5330,5366,5371,5372,4051,5010,5011,5012,5340,5343,5025,5046,5235,5251,5273,5364,4042,4044,4085,4086,4036,4029,4032,4026,4027,6507,4023,4024,6335,4249,4250,4621,4266,4267,4640,6292,6297,6884,6298,6301,6302,4048,4050,4052,6898,6899,6901,4045,4046,5228,4055,4059,5231,5236,4068,6894,6896] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.388,0.702]
select surf_pocket2, protein and id [2842,2874,237,239,240,271,3373,495,497,3191,3194,3159,3161,3372,3371,393,479,263,383,359,361,362,330,3163,242,244,246,3406,296,3173,245,270,2835,2837,2857,2841,2846,2847,2849,2850,2878,2882,2881,2931,3122,3124,2722,2687,2719,2728,367,364,366,368,369,3802,3803,331,2939,2940,2658,2685,499,500,681,2657,2659,680,2755,2762,2765,2752,2754] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.427,0.902]
select surf_pocket3, protein and id [6409,6411,6412,6413,6415,7528,7424,7033,7410,7440,7030,6976,6977,7524,6982,6831,6832,7439,6681,6687,6690,6694,6726,6725,6580,6581,6582,7527,7526,6653,6651,6654,7384,7385,7541,7519,7521,7522,7518,7418,7421,7536,7538,7413,7542,7549,7550,7551,7698,6683,7543,6682,6686,6375,6381,6430,6432,6433,6376,6379,6380,7003,7027,7026,6408,6417,6418,6406,6419,6420,6398,6385,6387,6389,6394,6395,6382,6377,7058,7062,6585,6579,6599,6587,6565,6994,6999] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.286,0.278,0.702]
select surf_pocket4, protein and id [6932,6934,6935,6936,6589,6590,6952,6604,4040,6591,6953,4039,6341,6442,6456,6441,6596,6586,6601,4286,4288,4289,4303,4311,4314,4290,4313,4076,4077,4078,4097,4038,4019,6339,6342,6295,6296,6297,4037,4044,4031,4049,6898,6901,6902,4053,4054,6895,6897] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.443,0.361,0.902]
select surf_pocket5, protein and id [6069,6071,6063,6070,2771,2773,2760,2889,2890,2913,3739,2887,2759,3736,2884,2891,2892,2784,2787,2807,2865,6000,6002,6009,6035,5998,6004,6003,5990,5976,5977,5979,3749,3754,3755,5971,3737,6027,2750,5915,5880,5912,5888,5914,5978,6093,2809,2810,2808,5859] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.404,0.278,0.702]
select surf_pocket6, protein and id [71,985,925,84,86,87,88,927,1803,1802,979,1260,1263,1293,1294,1291,1292,1766,986,1021,1022,1024,983,1027,1030,1053,1055,1016,1017,1020,131,134,127,111,117,118,120,1269,1271,1063,1065,1066,1421,58,949,950,914] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.592,0.361,0.902]
select surf_pocket7, protein and id [471,473,656,994,651,1038,410,1092,1080,1921,1922,1908,1913,1936,1944,664,1946,448,445,446,1940,1942] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.518,0.278,0.702]
select surf_pocket8, protein and id [2533,2550,4863,4871,4840,4841,4842,6110,6108,5161,4899,4862,4477,2486,2487,4898,2492,2515,2517,2534,2363,2364,6107,2362,2347,2468,2469,2345,4486,5177,5179,4475] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.741,0.361,0.902]
select surf_pocket9, protein and id [602,604,605,606,416,432,440,523,537,541,2628,483,486,516,517,518,414,2596,2602,2604,2607,439,532,619,527,610] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.635,0.278,0.702]
select surf_pocket10, protein and id [7551,7682,7687,7694,7698,7691,6682,7566,7567,6686,7543,6409,6411,6412,6413,6414,6415,6416,6656,6659] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.890,0.361,0.902]
select surf_pocket11, protein and id [4049,4060,6898,6901,6902,5306,6928,4343,4345,6921,6924,6895,6897,6297,6935,4312,4332,4314,4313] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.647]
select surf_pocket12, protein and id [3794,2952,2955,337,339,340,3134,3136,350,352,885,886,865,878,2980,2954,2970,822,3791,3145,832,864] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.757]
select surf_pocket13, protein and id [50,52,54,56,58,152,155,423,424,425,609,612,1063,1066,1083,1081,1029,1030,1062,618] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.533]
select surf_pocket14, protein and id [1198,1200,1219,1228,1251,1720,1719,1718,1741,1754,1794,1716,1815,1817] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.608]
select surf_pocket15, protein and id [7620,7841,7843,7844,7851,7854,7847,7849,7850,7876,8531,8512,8514,8516,8528,8508,7624,8497,7309,7282,8524,7308,7278,7281] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.416]
select surf_pocket16, protein and id [4419,4995,5223,5225,4417,4420,4987,4982,5242,5244,5465,5240,5259,4391,4392,4393,4394,4396,4398,4395,4397] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.459]
select surf_pocket17, protein and id [3524,3559,3550,3585,5984,2912,2895] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.298]
select surf_pocket18, protein and id [4887,4889,4891,4892,4898,2492,2501,2508,2510,2514,2515,2534,4465,4468,4476,4919,2326,4890] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.412,0.361]
select surf_pocket19, protein and id [6092,2802,2805,2801,2816,2818,6083,6057,6082,6084,6088,6089,6090,6091,2333,2339,2343,2345,2348,2352,6118,6106,6058] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.376,0.278]
select surf_pocket20, protein and id [7412,7050,7051,7054,7999,7090,7997,8011,8039,7953,8037,8046,7747,8052,8053,7411] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.561,0.361]
select surf_pocket21, protein and id [2127,2118,2119,708,2215,710,2156,1929,1937,1901,1931,2132,2133,2154,2155,2128,3854] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.494,0.278]
select surf_pocket22, protein and id [2962,3118,3332,3203,3204,2958,3136,3137,3205] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.710,0.361]
select surf_pocket23, protein and id [4076,4078,4095,4097,6339,6342,6344,4017,6340,4093,4096,4001,6341,6346] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.608,0.278]
select surf_pocket24, protein and id [458,459,2633,2635,446,1947,1963,1965,1967,1940,1942,1952,2613,2615,467,1944] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.859,0.361]
select surf_pocket25, protein and id [7786,7789,7792,7817,7756,7758,7795,8184,8186,7693,7706,7708,7722,8187] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.675,0.702,0.278]
select surf_pocket26, protein and id [2764,2927,2947,3785,2948,2974,3008,3768,3769] 
set surface_color,  pcol26, surf_pocket26 
   

deselect

orient
