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

load "data/AF-P08922-F1_clean.pdb", protein
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

load "data/AF-P08922-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [9103,10238,9108,9110,9113,9546,9117,9129,9550,9555,8358,8689,8360,8362,8693,8695,8697,8699,10235,8372,8703,9542,10253,9932,10254,9557,9922,9917,9918,10234,10236,10237,9909,9127,9135,9562,9559,9563,9134,9581,10295,8713,8714,8720,8375,8379,8718,8743,10306,8383,8386,10259,9933,9937,9939,10294,10298] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.392,0.702]
select surf_pocket2, protein and id [3614,3628,3269,3610,3613,3608,3609,3285,1921,1923,1907,1914,2259,3275,2925,2608,2268,2917,2921,2919,2609,2599,2298,3646,3649,3648,3658,2614,3626,2273,2269,2592,2593,2594,2903,2589,2591,2892,2902,3599,3588,3594,3592,3252,2254,1903,2247,1902,2581,3259,3244,3593,2924,2943,2944,2616,2928] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.431,0.902]
select surf_pocket3, protein and id [5357,6386,6387,5359,5685,6077,4715,6392,4709,5701,5695,5696,5699,5700,6086,6087,6088,6089,6093,6405,5024,5053,6459,4726,6460,4722,4724,4725,5365,5702,5718,6447,6448,6451,4737,6069,6073,5361,5350,5364,5011,5344] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.278,0.702]
select surf_pocket4, protein and id [14512,15059,15060,15064,15068,14509,14942,14943,14944,13946,14338,13899,13943,13891,13892,13893,13896,13897,13898,13945,13895,14076,14482,14465,14471,14930,14916,14516,13900,13930,14092,13901,13907,13908,14507,14496,14506,14488,14484,14487] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.431,0.361,0.902]
select surf_pocket5, protein and id [12087,12773,12086,12083,12084,12075,12076,3139,3141,2852,2846,2856,2855,3079,3080,3004,3006,2850,3044,3024,12098,12102] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.392,0.278,0.702]
select surf_pocket6, protein and id [8447,8410,8412,8413,8414,8390,8392,8393,8429,8517,8518,8734,11086,8406,10323,10325,8768,8932,8764] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.576,0.361,0.902]
select surf_pocket7, protein and id [8397,8398,10271,10278,10281,10264,10268,10273,10283,8423,8552,8422,8087,8094,7387,8198,8112,8113,8091,8199,8200,10266,8541,8546,8540,8279,8271,8283] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.506,0.278,0.702]
select surf_pocket8, protein and id [6310,6312,6299,6298,6307,6020,6230,6178,6181,6231,6022,6027,6194,6196,6018] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.722,0.361,0.902]
select surf_pocket9, protein and id [14155,14147,587,590,14140,14142,572,14144,570,574,575,581,458,393,394,413,415] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.616,0.278,0.702]
select surf_pocket10, protein and id [8910,9335,9007,9003,8921,8776,8925,10666,8777,8778,9351,10645,10646,9349,8798,8729,8757,9352,9374,9375,10651,9022,9023,8906,9012,9018] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.863,0.361,0.902]
select surf_pocket11, protein and id [3166,3084,3165,3167,2974,2975,2978,996,998,1342,3301,2935,2937,2994,3300] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.671]
select surf_pocket12, protein and id [6336,4547,6327,4551,4612,4614,4618,6368,6370,4545,4508,6347,4543,4541] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.792]
select surf_pocket13, protein and id [3807,4587,4588,4589,3805,3806,3787,3789,3790,3797,4574,4448,4450,4458,4459,4460,4411,4582,4431] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.561]
select surf_pocket14, protein and id [1973,1614,1618,1929,1611,1639,1640,1594,1944,1972,3636] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.647]
select surf_pocket15, protein and id [5386,5390,6793,6795,6790,6791,6797,6792,6796,5392,5522,7160,7162,7178,6798,6799,6801,6807] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.447]
select surf_pocket16, protein and id [8188,8136,8153,10199,10193,10196,10197,10189,10191,10217,8265,8207,10177] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.506]
select surf_pocket17, protein and id [3224,12505,12530,12506,12497,3033,12199,12511,2878] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.337]
select surf_pocket18, protein and id [5733,6117,6115,6119,6120,4397,5765,5711,5744,6098,6264,6100,6263] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.361]
select surf_pocket19, protein and id [5893,5902,5386,5390,5735,5896,5918,5920,5374,5375,5384,5389,5387,5719,5385,5716,5713,5911] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.337,0.278]
select surf_pocket20, protein and id [14,719,52,54,55,10,56,7,18,24,26,28,1,2,30,31,40,49] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.506,0.361]
select surf_pocket21, protein and id [6240,6248,6236,6306,3758,6234,6237,6251,6290,6293,6291,6130,6256,4379,3744,4384] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.447,0.278]
select surf_pocket22, protein and id [9238,9221,8891,9098,9219,9036,9039,8872,9062,8877,9065,9083] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.647,0.361]
select surf_pocket23, protein and id [12466,12468,12480,12475,12478,12479,12302,12303,12304,12518,12520,12445,12298,12308] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.561,0.278]
select surf_pocket24, protein and id [9146,8742,9580,9161,9162,10313] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.792,0.361]
select surf_pocket25, protein and id [4911,4555,4557,4630,4468,4469,4466,4467,4470,4471,4887,4903,4899,4900,4487,4648,4649] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.671,0.278]
select surf_pocket26, protein and id [3287,3301,2935,2936,2937,3293,2953,3300,2952,981,2948,2951,2954,2957,996,997,998,989,990] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.863,0.902,0.361]
select surf_pocket27, protein and id [12625,12628,12433,12452,12564,12581,12584,12431,12432,12435,12358,12355,12360] 
set surface_color,  pcol27, surf_pocket27 
   

deselect

orient
