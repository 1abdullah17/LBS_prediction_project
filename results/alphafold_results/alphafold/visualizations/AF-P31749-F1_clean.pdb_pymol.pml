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

load "data/AF-P31749-F1_clean.pdb", protein
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

load "data/AF-P31749-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [2188,1711,1709,1710,2084,2111,2112,2190,2191,1101,1088,1091,1093,1095,1100,1087,1141,3363,3372,3382,3377,3379,1081,1083,3417,2096,2195,2197,2199,2186,2192,2193,2198,2200,1263,1523,1524,2207,1446,1659,1660,1665,1262,1412,1138,1140,1246,1682,1686,1082] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.365,0.278,0.702]
select surf_pocket2, protein and id [1413,1104,1110,1118,1120,1263,1264,1280,1116,1122,2054,2195,2196,2197,2199,2198,1101,1100,1125,1141,1376,1375,1360,1383,1389,2214,2212,2213,2215,2220,2221,2222,1279,1384,2071] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.792,0.361,0.902]
select surf_pocket3, protein and id [2965,2290,2445,2450,2461,2974,2476,2479,2963,2276,2277,2279,2284,2462,2465,2468,2470,2471,2962,2964,3009,3010,2972,3011,3007,3015,3036,3038,2274,3014,3016] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.533]
select surf_pocket4, protein and id [1983,1947,1991,1925,1926,3203,3205,3197,3199,3201,3077,3157,3176,3158,1951,1491] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.361]
select surf_pocket5, protein and id [1727,1729,2086,2087,2088,2626,2356,2357,2359,2369,2371,2581,2582] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.533,0.278]
select surf_pocket6, protein and id [2046,2235,2315,2237,1379,1382,2222,1404,2332] 
set surface_color,  pcol6, surf_pocket6 
   

deselect

orient
