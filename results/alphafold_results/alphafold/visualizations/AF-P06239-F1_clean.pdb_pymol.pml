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

load "data/AF-P06239-F1_clean.pdb", protein
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

load "data/AF-P06239-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [2158,2180,2183,2154,2567,2568,2635,2186,2544,2552,2637,2543,2188,1787,1788,1789,2639,2640,2644,2645,2646,2668,1772,1784,2131,2132,1781,2116,1616,1617,1619,1620,1622,1670,1672,1673,2211,1631,1625,1627,1629,1630,1623,2004,2006,2005,2655,1927,2151,2153,2166] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.310,0.702]
select surf_pocket2, protein and id [190,1457,1496,188,191,1744,1718,1720,1745,1730,1743,1751,1752,1741,179,181,1742,180,1444,1452,1456,1433,1434,1440,1430,1431,1443,164,1435] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.498,0.361,0.902]
select surf_pocket3, protein and id [674,3508,3524,3532,3534,3535,3536,3543,3544,3548,3539,2337,2308,2332,3541,2307,2310,731,2320,984,691,693,694,696,996,990,695,1000,992,994,671,672,702,2344,2340,2342,2343,2361,3505,3500,2306,2312,2319] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.525,0.278,0.702]
select surf_pocket4, protein and id [1397,1417,628,624,560,559,561,627,140,141,128,131,251,252,569,493,130,1414,1416,564,127,716,1391,1392,619,621,584,582,576,606,614,715,1389,647] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.851,0.361,0.902]
select surf_pocket5, protein and id [1424,1425,1426,2146,1401,2169,2170,2587,2588,2589,2596,2160,2141,2624,2625,2579,2610,2608,2581,1976,1977] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.600]
select surf_pocket6, protein and id [682,2391,654,652,844,3511,3512,1980,1972,2392,2390,2389,3482] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.596]
select surf_pocket7, protein and id [2642,2646,2665,2643,2670,2671,2701,2511,2673,2734,2736,2733,2737,2738,2674,2515,2745,2552,2551,1652,1653,2721,2723,2728,2735,1641] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.325]
select surf_pocket8, protein and id [350,1465,1469,1485,1487,460,461,462,1474,1479,1486,196] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.478,0.361]
select surf_pocket9, protein and id [2239,3015,3027,3028,3032,3034,2530,2534,2776,2775,2540,2543,2202,2203,2978,3020] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.510,0.278]
select surf_pocket10, protein and id [1932,1933,1935,1936,1960,1962,1961,338,341,1447,1476,1478,1475,1472,335,2015,2018,346,353,2023,2031,1449,1461,1446] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.831,0.361]
select surf_pocket11, protein and id [1803,1838,2669,2692,2666,2689,2699,2701,1802,2665,1829,1646,1831] 
set surface_color,  pcol11, surf_pocket11 
   

deselect

orient
