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

load "data/AF-P00519-F1_clean.pdb", protein
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

load "data/AF-P00519-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [1933,2136,2121,1776,1653,1655,1598,1600,1604,1605,1608,1656,1780,2671,2670,2672,1610,1611,1612,1613,1640,1779,1907,2192,2562,2663,2589,2590,2195,2007,2570,2662,2665,2666,2576,2668,1616,1632,1625,2142,2157,2158,2162,2187,2189,2163,2186,2185,2137,1602,1764,1599,1603,2156] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.310,0.702]
select surf_pocket2, protein and id [664,847,1982,1984,673,2635,1992,1993,2358,2363,2365,2391,2394,2355,2361,1374,649,653,736,1371,1372,1349,1350,1351,1353,1354,735,737,1368,1373,1357,650,668,674,677,679,699,702,3594,694,846,848,672,2636,1380,1978,2412,2414,2381,2383,2386,2387,3560,3561,2382,2436,2384,3592,3590,2385] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.498,0.361,0.902]
select surf_pocket3, protein and id [2341,2344,3081,3083,3337,3084,3335,3689,3330,3331,3333,3344,3101,3323,2347,2352,3584,3368,3659,3369,3078,2377,2340,2342,3045,3080,2312,2314,2319,2320,2322,3688,3661,3664,3723,2346,2279,3103] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.525,0.278,0.702]
select surf_pocket4, protein and id [1405,193,199,201,203,188,1384,2171,1389,534,1380,1382,2172,2602,2605,2607,1746,1745,1724,1725,220,207,218,216,2139,2141,2147,2151,1402,1403,1421,2012,2014,1399,2652,1393,1397,1398,2636,2146,2164,2599,2166,1423,1427,2033,1751,1371] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.851,0.361,0.902]
select surf_pocket5, protein and id [1968,1410,364,372,375,382,374,395,396,379,1433,1446,1462,1432,1488,1492,1470,1473,1472,1474,1966,1967,1938,1951,1949,1463,1489,1487,1491,1950,1459] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.600]
select surf_pocket6, protein and id [683,685,686,824,3626,860,875,837,839,854,1014,1015,671,675,677,680,840,849,827,829,834,835,895,885,897,1011,1013,1008,989,3591,3627,3628,3630] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.596]
select surf_pocket7, protein and id [568,85,602,613,603,59,75,81,86,1176,1161,1163,1159,1164,583,585] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.325]
select surf_pocket8, protein and id [2712,2695,1870,2785,2801,2524,2525,2531,2694,2829,2807,2766,2775,2777,2820,2818] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.478,0.361]
select surf_pocket9, protein and id [1474,1497,1503,1701,1505,1510,1432,2030,1428,1434,222,1720,1718,1719,1732,1544,1426] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.510,0.278]
select surf_pocket10, protein and id [2301,2302,2303,2611,2230,2232,2625,2619,2627,2628,2262,2264,2182,2184,2181,2596,2601] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.831,0.361]
select surf_pocket11, protein and id [713,747,749,2317,955,970,973,750,957,961,967,2316,2318,2330,3665,714,3651,3663,715,717,3649] 
set surface_color,  pcol11, surf_pocket11 
   

deselect

orient
