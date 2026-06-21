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

load "data/AF-O15530-F1_clean.pdb", protein
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

load "data/AF-O15530-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [353,354,355,232,235,221,188,190,195,196,199,1291,1293,1294,1292,1296,1298,1299,520,234,337,774,777,1152,1168,1169,1193,1300,1451,1453,1454,1461,204,213,215,1313,1314,521,478,480,488,479,481,482,484,1321,1445,1449,487,460,461,828,830,180,792,178,183,184,186,177,800,783,1209,797,1208,823,801,820] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.278,0.702]
select surf_pocket2, protein and id [3417,3419,3443,2537,2517,3418,2515,2516,3441,2532,2826,3456,2831,2821,3454,3394,2839,3405,3407] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.576,0.361,0.902]
select surf_pocket3, protein and id [2301,604,605,615,1038,2317,2318,2275,2294,2302,2197,598,1037,1064,1057,1058,1060,1061,1062,595,1068] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.616,0.278,0.702]
select surf_pocket4, protein and id [1513,2073,1572,1574,1587,2071,2135,2136,1511,2057,2069,2072,2067,2079,1547,1394,1396,1546,1548,1578,1384,2124,2125,2127,1383,2114,2119,1584] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.792]
select surf_pocket5, protein and id [1880,1891,1894,1898,1877,1836,1747,3431,3433,3434,3442,3419,3443,3406,3411,1753,1750,1729,3249,3251,1874] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.447]
select surf_pocket6, protein and id [901,2386,3386,2851,2868,2870,2871,3362,3365,2911,3267,3266,3375,3382,2852,3289] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.361]
select surf_pocket7, protein and id [603,789,782,1276,600,602,2313,2316,1221,1223,1224,1242,1257,2314,621] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.447,0.278]
select surf_pocket8, protein and id [899,901,2386,3239,892,905,3386,2392,3267,3264,3266,3289] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.792,0.361]
select surf_pocket9, protein and id [977,981,2255,969,2258,2263,2256,2357,2358,973,2332,2337,2342,948,2308] 
set surface_color,  pcol9, surf_pocket9 
   

deselect

orient
