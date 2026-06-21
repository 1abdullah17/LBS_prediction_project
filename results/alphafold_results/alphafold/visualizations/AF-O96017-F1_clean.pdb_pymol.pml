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

load "data/AF-O96017-F1_clean.pdb", protein
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

load "data/AF-O96017-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [1195,1198,1201,1206,1208,2127,2152,2276,2128,2151,2107,2111,2299,2331,2295,2355,2356,2357,2109,2362,2368,2378,1180,1184,1185,1186,1188,1190,1192,1780,1782,1178,1194,1225,1227,1228,1349,1350,1731,2275,2267,2268,1179,1182,1332,1748,2140,2167,2166,1753,1754,1772,1750,1211,1212,1209,1203] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.329,0.278,0.702]
select surf_pocket2, protein and id [1799,1801,1802,1833,1822,1824,1827,1828,2648,2623,2624,2664,2666,2667,2411,2413,2401,1826,2386,2143,2650,2647] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.698,0.361,0.902]
select surf_pocket3, protein and id [1717,2298,2270,2277,2281,2324,2287,1525,2286,1348,1350,1731,2275,2269,2301,1716,1732,1616,2283,2285] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.639]
select surf_pocket4, protein and id [2550,3076,3078,3080,3082,3088,3096,3097,3057,3036,2495,2496,2498,2516,2517,2525,2526,2528,2043,2520] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.545]
select surf_pocket5, protein and id [1746,1312,3337,1763,1756,1741,3303,3305,3331,3306,3307,3308,2180,2179,2181,1736,1311] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.353,0.278]
select surf_pocket6, protein and id [1358,596,1221,1357,1065,1370,1372,1684,1685,1706,1709,593,594,1389,1677,1054,1055] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.729,0.361]
select surf_pocket7, protein and id [2298,2324,2320,2331,1206,1208,1366,2323,2332] 
set surface_color,  pcol7, surf_pocket7 
   

deselect

orient
