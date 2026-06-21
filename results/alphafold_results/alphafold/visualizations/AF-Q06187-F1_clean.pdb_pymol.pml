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

load "data/AF-Q06187-F1_clean.pdb", protein
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

load "data/AF-Q06187-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [3014,2969,2974,2975,3015,2964,3512,3515,3118,2961,2958,2957,2963,2962,3131,3128,3358,3495,3942,4013,4015,2994,3136,4022,4023,3489,3490,3134,3533,3540,3534,3537,3542,2971,3565,3516,2959,3943,3915,3929,3914] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.333,0.702]
select surf_pocket2, protein and id [3332,2781,2784,2785,2789,2790,3388,3498,4001,4002,3333,3496,3500,3504,3506,3091,3094,3089,3100,3103,3104,1507,2795,2796,3096,2778,2770,2776,2767,2768,2775,2764,3099,1480,1844,1847,1469,1485,3957,3954,3964,3965,3972,3334,3985] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.443,0.361,0.902]
select surf_pocket3, protein and id [3739,3334,3344,2035,3713,3715,2036,3769,3746,3747,1997,2735,2736,2737,1975,2762,3970,3986,3972,3767,3768,1989,1992,1990,1988,1996,2761,2768,2760,2756,3332,3333,3335] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.463,0.278,0.702]
select surf_pocket4, protein and id [1103,268,3277,3278,3266,3267,3273,271,1681,3298,3301,1093,24,1042,1672,23,270,9,10,1062,1063,1069,1041,1074,25,1073,269,1064,1101,1072,3269] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.741,0.361,0.902]
select surf_pocket5, protein and id [1803,1804,1806,2826,2829,1055,1084,1086,1385,1387,1134,1526,1527,1528,1031,1032,1067,1054,1056,1057,1058,1663,1665,1693,1386,1389,1650,1646,1717] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.694,0.278,0.702]
select surf_pocket6, protein and id [2801,2803,2838,3294,1699,2817,3302,1093,1695,2830] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.757]
select surf_pocket7, protein and id [3372,3375,3489,4046,3474,4045,3358,4013,4016,4017,4022,3286,3376] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.475]
select surf_pocket8, protein and id [1179,1330,1332,1141,1142,1156,1162,1339,883,885,1258,1261,1313,832,845,846,1173,1176] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.459]
select surf_pocket9, protein and id [4020,4125,4127,2994,2996,4019,4023,2981,2997,4042,3888,4050,2978,2975,3929,3917,3892,3928,3918] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.318,0.278]
select surf_pocket10, protein and id [1140,1372,1385,1387,1080,1082,1131,1134,510,506,2826,1055,1084,1086,1149,1150,1151,1056,1058,1054,515] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.561,0.361]
select surf_pocket11, protein and id [2185,4956,4959,2173,2174,2180,2182,2184,4960,1989,1991,1990,1993,1987,2000,1994,4990,4991,2169,2009] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.553,0.278]
select surf_pocket12, protein and id [3313,3315,3318,3316,3311,3312,3317,3801,1676,1682,1683,1685,1688,1670,3299,1827,1831] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.859,0.361]
select surf_pocket13, protein and id [109,858,875,458,447,451,313,316,444,110] 
set surface_color,  pcol13, surf_pocket13 
   

deselect

orient
