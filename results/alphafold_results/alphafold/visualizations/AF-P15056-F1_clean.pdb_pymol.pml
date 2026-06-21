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

load "data/AF-P15056-F1_clean.pdb", protein
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

load "data/AF-P15056-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [2821,2712,2723,2741,2822,2826,2828,1931,2827,2739,2740,2742,1830,1833,1785,1787,1788,1914,2298,2320,2321,2312,2317,2293,1832,2315,2316,1783,1779,2818,2837,2173,2176,2177,2174,2292,2737,2342,2738,2335,2337,2345,2347,2711,1929,1930,2850,1926,2098,2097,2279,2278] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.329,0.278,0.702]
select surf_pocket2, protein and id [1202,1203,1394,1396,1398,2624,1206,1208,1472,1473,1537,1535,1536,1543,1513,1514,1491,2614,2639,3028] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.698,0.361,0.902]
select surf_pocket3, protein and id [2597,2143,2146,2152,2138,1560,622,623,621,2589,2623,1179,830,1165,1559,1561,1174,812,822,2599,2632,2130,1162,2568,1163] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.639]
select surf_pocket4, protein and id [1436,1437,1189,1544,1449,1451,1590,1592,1594,1443,1558,1542,1571,1586,1568,1574] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.545]
select surf_pocket5, protein and id [2807,2757,2760,2768,2769,2791,2792,632,2300,2302,2301,2808,2181,2185,2153,620,627,638,2148,2141,640,630] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.353,0.278]
select surf_pocket6, protein and id [1796,2723,2828,1814,2684,2700,2701,2858,1804,1812] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.729,0.361]
select surf_pocket7, protein and id [1806,1811,2008,1962,1809,1943,1946,1947,2039,2043,2009,2072,2851,2848] 
set surface_color,  pcol7, surf_pocket7 
   

deselect

orient
