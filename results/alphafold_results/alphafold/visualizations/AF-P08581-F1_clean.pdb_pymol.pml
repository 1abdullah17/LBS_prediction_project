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

load "data/AF-P08581-F1_clean.pdb", protein
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

load "data/AF-P08581-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [7486,7487,7610,7626,7627,7983,8420,8499,8501,8502,8504,7863,7760,8392,8508,7442,7437,8008,8031,8032,8007,8419,8038,8418,7436,7438,7440,8042,8004] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.365,0.702]
select surf_pocket2, protein and id [516,520,524,525,889,496,503,505,613,925,927,910,924,926,1074,1076,721,491,866,880,882,886,722,867,594,600,602,607,608,596,610,593,858] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.361,0.902]
select surf_pocket3, protein and id [239,234,238,756,3587,243,3586,1243,1809,1812,1241,1242,1800,1804,1813,1248,766,748,752,1816,244,3592,3600,3576,3588,2331,1811,1817,1821,2339,2341,2343,2344,2347,2348] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.365,0.278,0.702]
select surf_pocket4, protein and id [2154,2772,3022,3024,1868,2133,2971,2972,2973,2977,2978,2979,2980,3021,2059,2062,2755,2134] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.576,0.361,0.902]
select surf_pocket5, protein and id [3989,4144,4329,4148,4162,3863,3865,4323,4325,4327,4341,3868,3973,3975,3977,3857,58,60,4333] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.533,0.278,0.702]
select surf_pocket6, protein and id [4595,4597,3867,3876,3877,3880,3881,3885,3878,3908,4418,3906,4588,4590,4591,4592,65,64,66,4607] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.792,0.361,0.902]
select surf_pocket7, protein and id [1224,1223,874,1213,1216,902,892,897,1428,1451,891,893,1430,1396,1406,1407,1212,1470,1204,1469] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.702]
select surf_pocket8, protein and id [3989,4144,4147,4150,4148,4162,3971,3973,3975,3977,3871,4134,4139,4110,4341] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.792]
select surf_pocket9, protein and id [3466,3411,3917,3919,3294,3413,3414,3462,3453,3457,3459,3461,3468,3907,3915,2389,2390,3425,3427,3295] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.533]
select surf_pocket10, protein and id [3467,3468,3906,3850,3852,3905,3907,3910,3848,3909,3881,3885,3908,65,79,83,3466,3479,3465,3473] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.576]
select surf_pocket11, protein and id [5056,5058,4711,4709,4854,4856,4874,5225,5075,5076,5074,5077,4221,4231,4229,4220,4222] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.365]
select surf_pocket12, protein and id [8078,8080,8130,8462,8463,8464,8465,8028,8029,8133,8121,8125,8127,8134,8148,8136,8138,8142] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.361]
select surf_pocket13, protein and id [5293,7668,5296,5308,7938,7939,5309,5292,7198,7942,5160,5161,7933,7945,7947,7638,5159] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.365,0.278]
select surf_pocket14, protein and id [1243,1242,1782,225,237,238,1222,1227,745,1230,1234,1788,743,227,746] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.576,0.361]
select surf_pocket15, protein and id [8297,9250,8710,9227,9231,9236,9233,8677,8678,8679,8680,8687,8688] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.533,0.278]
select surf_pocket16, protein and id [4497,4505,4508,4516,4518,4502,4666,4660,4664,4696,4698,4699,4700,4701,4702,4728,4729,4864,4680,4866] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.792,0.361]
select surf_pocket17, protein and id [1455,1095,1105,580,582,1442,1443,1454,1094,1100,1109,1112,1113,1161,1167,1171,1148] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.702,0.278]
select surf_pocket18, protein and id [506,507,721,504,732,867,866,734,737] 
set surface_color,  pcol18, surf_pocket18 
   

deselect

orient
