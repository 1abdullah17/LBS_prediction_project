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

load "data/AF-P52333-F1_clean.pdb", protein
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

load "data/AF-P52333-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [3831,3815,3814,3806,3809,3812,3832,4764,3799,4431,4432,4924,3822,3824,4904,4349,4354,4019,4020,4002,4003,4787,3849,4379,4400,4401,4386,4921,4248,4359,3787,3786,3789,3790,4373,4771,4772,4773,4748,4758,4759,4760,4788,4404,4407,4734,4928,4735,5025,5028,4761,5023] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.396,0.702]
select surf_pocket2, protein and id [640,643,1744,1786,1789,1790,1793,1788,1800,1740,1771,1774,1742,635,637,614,1806,1808,639,606,611,598,931,967,411,437,998,435,408,409,906,940,941,942,1004,1005,1007,1009,1010,995,965,963,999,601,1794,1796,1795,1798] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.447,0.902]
select surf_pocket3, protein and id [6936,6944,6965,6966,7040,6950,6935,7042,5996,6121,6524,6525,5997,7043,7045,6138,7047,7049,6136,6137,6913,5951,5958,5959,5946,5947,5949,5948,5994,5984,5973,5977,5981,5979,5968,6551,5940,5943,5944,6542,6544,5941,6569,6570,6572,6604,6575,6577,6559,6381,6547] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.298,0.702]
select surf_pocket4, protein and id [6554,6979,6980,6981,6982,6535,6539,2308,2309,2312,6107,4201,6101,6103,7006,6986,6988,6389,6534,6527,6529,6413,6410,7028,4656,4658,2319,4209,4196,4200,4203,4691,4693,4694,2304,2303,2301,2306,2302,4681,4685,7010,6353,6354,6355,7011,2104,2107,6390] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.404,0.361,0.902]
select surf_pocket5, protein and id [334,356,364,366,337,368,463,458,460,462,332,338,339,306,2688,2634,2656,212,2648,2745,2746,2929,342,2552,2481,2513,2514,2517,2636,2456,2455,341,343,2523,2525,3552,3551,3554,630,631,2689,2659,2736,2737,233,235,336,2655,2510,2551,2509,2511] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.365,0.278,0.702]
select surf_pocket6, protein and id [2271,2336,2317,2129,2286,4623,4625,4627,4593,4592,2161,4854,2353,2337,4570,4568,4599,4601,4591,4569,2144,2160,5761,5782,2269,2272,4233,4870,4872,4234,4224] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.533,0.361,0.902]
select surf_pocket7, protein and id [2558,5679,5684,5710,5737,5745,5744,3505,3506,2557,2562,3486,5649,5656,5665,5673,5675,5595,5596,5597,2563,2564,2566,2573,2546,2556,2547,2570,2668,2695,2696,2581,2583,2671,2608,2694] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.467,0.278,0.702]
select surf_pocket8, protein and id [4178,3662,4177,6559,4170,4262,4256,4257,3678,3681,4270,4267,4269,4273,5920,6643,3672,3675,3677,3679,3680,4277,5944,6546,6014,5941,5932,5938,5927,5925,4193] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.663,0.361,0.902]
select surf_pocket9, protein and id [6101,6036,6057,6058,6098,4224,6095,2309,6107,6116,6038,6040,2328,5843,5844,5845,6043,6413,6409,6410,5814,5816,5817,6407,5831,4217,4220,4221,6091,2321,2326,2315,2317,2316] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.565,0.278,0.702]
select surf_pocket10, protein and id [6348,6342,6343,2052,2060,2067,2053,1974,2297,2063,2329,6353,6354,6355,2076,2107,2106,2068,2069,2073,2054,6338,6341,2048,6334,6336,6333,6328,6391,6394,6395,6397,6407,5825,6383,6385,6387,6388,6389,6413,6409,2328,6346,6314] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.792,0.361,0.902]
select surf_pocket11, protein and id [4711,4712,4932,4934,4936,5148,5144,4718,4918,4140,4141,4947,4954,3602,4159,4161,4164,4713,3609,3611,3613,4187,4706,4142,4143,4145,4136,4923,4139,4158] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.667,0.278,0.702]
select surf_pocket12, protein and id [6236,6204,6232,6234,6201,6203,6205,6207,5978,7072,7209,7212,7205,7195,7200,7196,7179,5973,5976,5974,6155,7085,7163,6905,7086,6231,6237,7064,7066,6276,6238] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.878]
select surf_pocket13, protein and id [5021,5023,5024,5047,5053,5051,5052,5056,5057,5267,5305,5025,4463,4460,5014,5337,5339,5309,4748,4746,4754,4757,4760,4420,4424,4735,3814,5311] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.631]
select surf_pocket14, protein and id [763,814,759,761,766,3576,740,3588,3590,3592,3593,738,703,739,777,779,786,791,792,793,794,3582,8069,8107,8109,8070] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.749]
select surf_pocket15, protein and id [819,821,807,3543,2758,2760,2761,2787,2777,327,2943,2954,314,325,315,2944,756] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.533]
select surf_pocket16, protein and id [2344,1852,1862,1863,1861,1974,2053,2063,2329,2065,2334,2339,2346,1971,1972,2048,1973,1976,1977,5818,5816,2328] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.620]
select surf_pocket17, protein and id [63,106,395,966,959,1000,1003,393,928,925,957,62,57,59,583,44,596,601,56,577,578,42,43,46,50,52,988,986] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.431]
select surf_pocket18, protein and id [3355,3356,3357,3360,3361,2896,3465,2706,2717,2719,3295,3298,3299,3325,3328,3324,3330,3331,5659,3480] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.490]
select surf_pocket19, protein and id [4275,3974,3978,3980,3984,3983,3985,3988,4367,6086,6084,6085,4368,4356,4358,4360,4362,4364,6031,4257,5919,5896,5899,5906] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.329]
select surf_pocket20, protein and id [7556,7569,7554,7266,7522,7537,6592,7519,7521,7520,7478,7565,7265,6596,6938] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.361,0.361]
select surf_pocket21, protein and id [2733,2736,2737,2753,2689,3554,3519,3524,2817,3521,3522,3493,2559,2560,2562,3520,2553,3539,3537,3529] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.329,0.278]
select surf_pocket22, protein and id [456,2402,2404,445,478,477,424,2400,2426,447] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.490,0.361]
select surf_pocket23, protein and id [2588,2596,2598,2149,5787,2135,2140,2142,2145,5789,5757,5756,5758,5762,5763,5765,5731,5732,2148,2150,2151,2152,2397,2615,2616,2618,2617,2601,2604] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.431,0.278]
select surf_pocket24, protein and id [3657,3660,6639,4179,6559,6643,6564,6566,6563,4185,6567,6614,6565] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.620,0.361]
select surf_pocket25, protein and id [4187,3610,6996,6680,4188,3660,4157,6653,3641] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.533,0.278]
select surf_pocket26, protein and id [7072,7208,7211,7216,7226,7196,7179,6905,7180,7276,7224,7240,7695,7335,7310,7274] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.749,0.361]
select surf_pocket27, protein and id [7933,7934,7414,7792,7935,7936,7797,7802,7596] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.631,0.278]
select surf_pocket28, protein and id [1323,1324,1462,1465,1510,1414,1415,1441,1443,1469,1507] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.902,0.878,0.361]
select surf_pocket29, protein and id [7358,7360,7377,7963,7965,7966,7146,7148,7149,7151,7342,7343,7294] 
set surface_color,  pcol29, surf_pocket29 
set_color pcol30 = [0.667,0.702,0.278]
select surf_pocket30, protein and id [644,673,639,641,648,625,633,842] 
set surface_color,  pcol30, surf_pocket30 
   

deselect

orient
