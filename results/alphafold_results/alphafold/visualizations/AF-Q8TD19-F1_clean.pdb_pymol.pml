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

load "data/AF-Q8TD19-F1_clean.pdb", protein
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

load "data/AF-Q8TD19-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [2593,2585,5161,5152,3449,3856,3857,3858,4350,4345,4355,4365,3435,2577,2579,3024,3025,3441,4747,5166,5174,4357,4363,2598,4730,5146,4739,5153,3044,3863,3039,3454] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.341,0.702]
select surf_pocket2, protein and id [699,690,536,698,1125,1126,1129,1128,1130,256,258,259,672,149,723,720,695,678,277,279,280,281,282,1208,1205,1206,1210,1212,535,432,82,692,694,97,93,96,92,94,98,99,712,713,714,715,717,752,788,705] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.424,0.361,0.902]
select surf_pocket3, protein and id [1158,2225,3124,3126,2245,2247,708,1142,1144,1145,1148,1146,1174,2686,2688,2997,3118,707,2689,3114,3116,3117,706,2254,2252,2246,686,681,687,242,2275,228,238,235,2669,241,243,244,689,704] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.435,0.278,0.702]
select surf_pocket4, protein and id [5244,2474,2407,2408,2410,2411,2413,2472,2392,2394,2393,2434,4859,5228,5243,5080,5082,5083,5084,5114,5115,2379,5112,5097] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.698,0.361,0.902]
select surf_pocket5, protein and id [504,501,990,490,965,986,988,1017,1015,2033,984,508,510,512,936,966,2175,2177,2178,2179,2058,2139,2140,2143,2202,2204,520,934,2176,2173,2212] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.651,0.278,0.702]
select surf_pocket6, protein and id [1226,1227,1234,1228,369,399,407,433,400,405,398,404,1071,1210,1211,1212,125,129,340,339,281,282,298,297] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.824]
select surf_pocket7, protein and id [2018,1001,1002,1006,1477,1289,1481,1491,1492,2020,1003,1009,1010,2034,1483,1286,1287,1471,1472] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.533]
select surf_pocket8, protein and id [156,144,157,290,291,2353,5300,632,143,145,2340,273,177,60,2350,2352,2354,5299,631] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.545]
select surf_pocket9, protein and id [2744,2756,2760,2761,2455,2736,2897,2501,2896,2749,2420,2479,2751] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.318]
select surf_pocket10, protein and id [2607,2608,2610,2770,2771,5191,5193,5195,4378,5183,5189,4761,4765,5180,5203,5206] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.451,0.361]
select surf_pocket11, protein and id [22,25,56,205,5269,5119,4840,5263,5266,188,5121,5123,5124,5248,5253,5122,5125] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.459,0.278]
select surf_pocket12, protein and id [664,557,663,665,569,572,196,199,2329,559,560,2271,558,2280,46,48,214,213,215] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.729,0.361]
select surf_pocket13, protein and id [491,493,2265,2269,495,540,2232,2234,2236,497,499,511,525,682,1192,2257,2261,562,2255,2253,2256] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.675,0.278]
select surf_pocket14, protein and id [2979,2936,2937,2696,2698,2697,82,2678,2679,2680,2529,2538,2549] 
set surface_color,  pcol14, surf_pocket14 
   

deselect

orient
