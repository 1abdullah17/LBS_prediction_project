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

load "data/AF-P12931-F1_clean.pdb", protein
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

load "data/AF-P12931-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [2571,1554,2501,2486,1608,1555,1550,1551,1552,2578,2579,2580,1720,2602,1721,1722,1772,1582,1737,1586,1717,1558,1609,1719,1592,1593,1590,2092,2112,2115,2502,2066,2051,2065,1704,1944,2573,2599,2604,2607,2605,2635,1865,1958,1943,2589,2478,2118,2144,1557,1560,1562,1569,1570,1573,1578,1588,2485,2449,2468,2477,2727,2679] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.310,0.702]
select surf_pocket2, protein and id [1375,2094,1351,2519,2529,1374,1373,1349,1356,1345,1346,1350,1376,2075,2080,2514,2518,2526,2511,1917,1930,2559,2560,2542,2544,1916,1921,629,630,631,666,2545,665] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.498,0.361,0.902]
select surf_pocket3, protein and id [92,103,106,107,1365,1369,1373,1377,454,457,98,100,1382,1355,1358,1360,1356,87,1375,1687,1681,1673,125,128,1972,1974,1383,1404] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.525,0.278,0.702]
select surf_pocket4, protein and id [2674,2676,2678,2741,2742,2861,2862,2813,2863,2864,2840,2842,2839,2683,2685,2687,2435,2689,2692,2818,2608,2610,2632,2434,2437,2607,2606] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.851,0.361,0.902]
select surf_pocket5, protein and id [69,70,529,532,509,555,557,523,592,594,564,565,567,568,569,570,571,572,573,590,442,1366,1345,1350,1367,1344,80,82,1364] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.600]
select surf_pocket6, protein and id [681,682,612,642,645,646,647,3543,616,618,619,620,623,3502,2277,2279,2281,2282,953,937,942,949,952,2249,2276,3494,3534,3530,3529,3531,3538] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.596]
select surf_pocket7, protein and id [3472,3473,1918,1919,1929,1909,1912,1920,2301,2303,2327,2329,2299,2302,2382,2385,600,601,624,632,616,618,619,3504,3505,3506] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.325]
select surf_pocket8, protein and id [3006,3026,3023,3024,3035,3062,3063,3210,3211,3200,3202,3203,3039,3005,3231,3221,3219,3223,2975,2972,2974] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.478,0.361]
select surf_pocket9, protein and id [1953,1392,1899,302,1902,1951,1870,1871,1873,1874,1967,1394,1409,1431,1435,1429,1959,1386,1388,1393] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.510,0.278]
select surf_pocket10, protein and id [2135,2136,2172,1578,2983,2779,2780,2743,2744,2468,2477,2746,2727,3036,3037,3024,3032] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.831,0.361]
select surf_pocket11, protein and id [3056,3192,3340,3337,3339,2922,2952] 
set surface_color,  pcol11, surf_pocket11 
   

deselect

orient
